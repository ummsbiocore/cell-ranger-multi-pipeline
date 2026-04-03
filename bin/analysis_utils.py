#!/usr/bin/env python

import scanpy as sc
import scvelo as scv
import cellrank as cr
import pandas as pd
import numpy as np
import os
from scipy import sparse
from scipy.stats import ks_2samp, mannwhitneyu

def get_optimal_threads(n_cells):
    """
    Determine optimal number of threads based on cell count.
    Ranges:
    < 5000: 1
    5000 - 10000: 4
    10000 - 20000: 8
    20000 - 30000: 12 (Interpolated)
    30000 - 40000: 16
    40000 - 50000: 32
    >= 50000: 48
    """
    if n_cells < 5000:
        threads = 1
    elif n_cells < 10000:
        threads = 4
    elif n_cells < 20000:
        threads = 8
    elif n_cells < 30000:
        threads = 12
    elif n_cells < 40000:
        threads = 16
    elif n_cells < 50000:
        threads = 32
    else:
        threads = 48
    
    # Cap at system limit
    try:
        cpu_count = os.cpu_count() or 1
        return min(threads, cpu_count)
    except Exception:
        return threads

def run_velocity_pipeline(subset, n_jobs=None):
    """Recompute velocity information on a subset copy."""
    # Work on a copy so we do not mutate the cached base AnnData
    data = subset.copy()

    # Guard against tiny subsets where neighborhood graphs cannot be computed
    if data.n_obs <= 2:
        raise ValueError("Not enough cells to recompute velocities for this condition.")

    # Determine parameters consistent with process_adata.py
    n_neighbors = max(5, min(30, data.n_obs - 1))
    n_pcs = 30
    if 'X_pca' in data.obsm:
        n_pcs = min(n_pcs, data.obsm['X_pca'].shape[1])

    # Recompute moments and velocity graph for the subset
    # Drop any pre-existing neighbor graph inherited from the full dataset
    data.uns.pop('neighbors', None)
    data.obsp.pop('distances', None)
    data.obsp.pop('connectivities', None)

    sc.pp.neighbors(data, n_neighbors=n_neighbors, n_pcs=n_pcs)
    scv.pp.moments(data, n_pcs=n_pcs, n_neighbors=n_neighbors)
    scv.tl.velocity(data)
    
    # Determine threads for this subset if not provided
    if n_jobs is None:
         n_jobs = get_optimal_threads(data.n_obs)
         
    scv.tl.velocity_graph(data, n_jobs=n_jobs)
    
    # Compute confidence if possible
    try:
        scv.tl.velocity_confidence(data)
    except Exception:
        pass

    return data

def compute_pair_probabilities(T, cluster_indices, filter_set):
    """Derive per-cell transition probabilities grouped by cluster pairs."""
    pair_probs = {}
    ordered_clusters = sorted(cluster_indices.keys())

    for src in ordered_clusters:
        if filter_set and src not in filter_set:
            continue

        src_idxs = cluster_indices[src]
        # No minimum cell filtering - include all clusters like CellRank does
        
        src_slice = T[src_idxs, :]
        for tgt in ordered_clusters:
            if filter_set and tgt not in filter_set:
                continue

            tgt_idxs = cluster_indices[tgt]
            if len(tgt_idxs) == 0:
                continue

            sub_matrix = src_slice[:, tgt_idxs]
            probs = np.asarray(sub_matrix.sum(axis=1)).flatten()
            pair_probs[(src, tgt)] = probs

    return pair_probs

def perform_transition_analysis(
    data,
    group_col,
    k_list,
    filter_set,
    n_macrostates,
    condition_col=None,
    compute_condition_cmp=False,
    n_jobs=1
):
    """Run transition, PAGA, and CellRank analyses on the provided AnnData."""
    # Enforce determinism
    np.random.seed(0)
    scv.settings.seed = 0

    a = data.copy()

    # Parse k list as integers with fallback
    try:
        k_values = [int(k) for k in k_list]
    except Exception:
        k_values = [1, 2, 3]

    result = {
        'group_col': group_col,
        'condition_col': condition_col,
        'k_list': k_values,
        'cluster_filter': list(filter_set) if filter_set else []
    }

    # Ensure neighbors exist for downstream PAGA / CellRank computations
    if 'neighbors' not in a.uns:
        pcs_available = a.obsm['X_pca'].shape[1] if 'X_pca' in a.obsm else 30
        sc.pp.neighbors(a, n_neighbors=min(30, max(5, a.n_obs - 1)), n_pcs=min(30, pcs_available))

    # PAGA trajectory
    try:
        # sc.tl.paga(a, groups=group_col)
        # Ensure connectivity keys exist
        if 'distances' in a.obsp:
            a.uns['neighbors'] = a.uns.get('neighbors', {})
            a.uns['neighbors']['distances'] = a.obsp['distances']
        if 'connectivities' in a.obsp:
            a.uns['neighbors'] = a.uns.get('neighbors', {})
            a.uns['neighbors']['connectivities'] = a.obsp['connectivities']
            
        scv.tl.paga(a, groups=group_col)
        result['paga_result'] = {'computed': True}
    except Exception as e:
        result['paga_result'] = {'computed': False, 'error': str(e)}

    # Velocity graph / transition matrix
    if 'velocity_graph' not in a.uns:
        scv.tl.velocity_graph(a, n_jobs=n_jobs)

    T = scv.utils.get_transition_matrix(a)
    if not sparse.issparse(T):
        T = sparse.csr_matrix(T)

    clusters = a.obs[group_col].astype(str)
    cluster_labels = clusters.values
    unique_clusters = np.unique(cluster_labels)
    cluster_indices = {str(c): np.where(cluster_labels == c)[0] for c in unique_clusters}

    result['cluster_labels'] = cluster_labels
    result['unique_clusters'] = unique_clusters
    result['cluster_indices'] = cluster_indices
    result['T'] = T

    # Compute cluster-level transitions using CellRank coarse graining when possible
    cluster_matrix_by_k = {}
    coarse_success = False
    try:
        vk_temp = cr.kernels.VelocityKernel(a)
        vk_temp.compute_transition_matrix(n_jobs=n_jobs)
        gpcca_temp = cr.estimators.GPCCA(vk_temp)
        n_states_temp = min(len(unique_clusters), n_macrostates)
        gpcca_temp.fit(cluster_key=group_col, n_states=n_states_temp) # n_jobs removed as per previous fix
        coarse_T = gpcca_temp.coarse_T
        coarse_success = True

        for k in k_values:
            if k == 1:
                mat = coarse_T.copy()
            else:
                mat = coarse_T.values.copy()
                for _ in range(k - 1):
                    mat = mat @ coarse_T.values
                mat = pd.DataFrame(mat, index=coarse_T.index, columns=coarse_T.columns)

            cluster_matrix_by_k[int(k)] = pd.DataFrame(mat.values, index=mat.index.astype(str), columns=mat.columns.astype(str))

    except Exception as coarse_exc:
        # Fallback to manual aggregation on full transition matrix
        print(f"[Transition Analysis] CellRank coarse-grained transitions failed: {coarse_exc}\nFalling back to direct aggregation.")
        cluster_list = sorted(cluster_indices.keys())
        for k in k_values:
            T_k = T.copy()
            for _ in range(int(k) - 1):
                T_k = T_k.dot(T)

            rows = []
            for src in cluster_list:
                src_idxs = cluster_indices[src]
                # No minimum cell filtering - include all clusters

                row_vals = []
                src_slice = T_k[src_idxs, :]
                for tgt in cluster_list:
                    tgt_idxs = cluster_indices[tgt]
                    if len(tgt_idxs) == 0:
                        row_vals.append(0.0)
                        continue
                    sub_matrix = src_slice[:, tgt_idxs]
                    probs = np.asarray(sub_matrix.sum(axis=1)).flatten()
                    row_vals.append(float(probs.mean()) if len(probs) else 0.0)

                rows.append(row_vals)

            cluster_matrix_by_k[int(k)] = pd.DataFrame(rows, index=cluster_list, columns=cluster_list)

    # Apply optional cluster filter
    if filter_set:
        for k, mat in list(cluster_matrix_by_k.items()):
            if mat.empty:
                continue
            rows = [idx for idx in mat.index if idx in filter_set]
            cols = [col for col in mat.columns if col in filter_set]
            cluster_matrix_by_k[k] = mat.loc[rows, cols] if rows and cols else pd.DataFrame()

    if coarse_success:
        print("[Transition Analysis] Coarse-grained cluster transitions computed via CellRank GPCCA.")

    result['cluster_matrix_by_k'] = cluster_matrix_by_k
    result['pair_probabilities'] = compute_pair_probabilities(T, cluster_indices, filter_set)

    # Terminal states and fate probabilities via CellRank
    fate_result = {'computed': False}
    df_condition_comp = None
    df_lineage_paths = pd.DataFrame()

    try:
        scv.tl.terminal_states(a)

        vk = cr.kernels.VelocityKernel(a)
        vk.compute_transition_matrix(n_jobs=n_jobs)
        estimator = cr.estimators.GPCCA(vk)
        n_states = min(len(unique_clusters), n_macrostates)
        estimator.fit(cluster_key=group_col, n_states=n_states)
        estimator.predict_terminal_states()
        
        try:
            # Try default computation first
            estimator.compute_fate_probabilities(n_jobs=n_jobs)
        except Exception as e:
            print(f"[Transition Analysis] Default fate computation encountered numerical instability. Switching to robust solver (gmres)...")
        
        # Check if results are present; if not, retry with robust settings
        if 'lineages_fwd' not in a.obsm:
            try:
                estimator.compute_fate_probabilities(n_jobs=n_jobs, solver='gmres', tol=1e-6, preconditioner='ilu')
                print("[Transition Analysis] Robust fate computation succeeded.")
            except Exception as e:
                print(f"[Transition Analysis] Robust fate computation also failed: {e}")

        try:
            estimator.compute_absorption_times()
        except Exception:
            pass
        
        fate_probs = a.obsm['lineages_fwd']
        lineage_names = fate_probs.names if hasattr(fate_probs, 'names') else []

        # Condition comparison
        if compute_condition_cmp and condition_col and len(lineage_names) > 0:
            conditions = a.obs[condition_col].unique()
            if len(conditions) >= 2:
                comparison_rows = []
                for lineage_name in lineage_names:
                    lineage_idx = list(lineage_names).index(lineage_name)
                    fate_vals = fate_probs[:, lineage_idx]
                    if hasattr(fate_vals, 'X'):
                        fate_vals = fate_vals.X.flatten()
                    else:
                        fate_vals = np.asarray(fate_vals).flatten()
                    for i in range(len(conditions)):
                        for j in range(i + 1, len(conditions)):
                            cond1, cond2 = conditions[i], conditions[j]
                            mask1 = a.obs[condition_col] == cond1
                            mask2 = a.obs[condition_col] == cond2
                            vals1 = fate_vals[mask1]
                            vals2 = fate_vals[mask2]
                            if len(vals1) == 0 or len(vals2) == 0:
                                continue
                            ks_stat, ks_p = ks_2samp(vals1, vals2)
                            mw_stat, mw_p = mannwhitneyu(vals1, vals2, alternative='two-sided')
                            comparison_rows.append({
                                'terminal_state': str(lineage_name),
                                'condition_1': cond1,
                                'condition_2': cond2,
                                'fate_prob_mean_cond1': float(vals1.mean()),
                                'fate_prob_mean_cond2': float(vals2.mean()),
                                'fate_prob_std_cond1': float(vals1.std()),
                                'fate_prob_std_cond2': float(vals2.std()),
                                'ks_statistic': float(ks_stat),
                                'ks_pvalue': float(ks_p),
                                'mannwhitneyu_statistic': float(mw_stat),
                                'mannwhitneyu_pvalue': float(mw_p)
                            })

                if comparison_rows:
                    df_condition_comp = pd.DataFrame(comparison_rows)
                    from statsmodels.stats.multitest import multipletests
                    if len(df_condition_comp) > 0:
                        _, ks_p_fdr, _, _ = multipletests(df_condition_comp['ks_pvalue'], method='fdr_bh')
                        _, mw_p_fdr, _, _ = multipletests(df_condition_comp['mannwhitneyu_pvalue'], method='fdr_bh')
                        df_condition_comp['ks_pvalue_fdr'] = ks_p_fdr
                        df_condition_comp['mannwhitneyu_pvalue_fdr'] = mw_p_fdr

        # Lineage paths summary
        path_rows = []
        for terminal_state in lineage_names:
            lineage_idx = list(lineage_names).index(terminal_state)
            fate_vals = fate_probs[:, lineage_idx]
            if hasattr(fate_vals, 'X'):
                fate_vals = fate_vals.X.flatten()
            else:
                fate_vals = np.asarray(fate_vals).flatten()

            for cluster in unique_clusters:
                cluster_mask = cluster_labels == cluster
                if cluster_mask.sum() == 0:
                    continue
                mean_fate = float(fate_vals[cluster_mask].mean())
                std_fate = float(fate_vals[cluster_mask].std())
                if mean_fate > 0.01:
                    path_rows.append({
                        'source_cluster': str(cluster),
                        'terminal_state': str(terminal_state),
                        'fate_probability_mean': mean_fate,
                        'fate_probability_std': std_fate
                    })

        if path_rows:
            df_lineage_paths = pd.DataFrame(path_rows).sort_values('fate_probability_mean', ascending=False)

        fate_result = {
            'computed': fate_probs is not None,
            'terminal_states': list(lineage_names) if fate_probs is not None else [],
            'n_macrostates': len(estimator.macrostates.cat.categories) if hasattr(estimator, 'macrostates') else 0
        }

    except Exception as fate_exc:
        fate_result = {'computed': False, 'error': str(fate_exc)}
        print(f"[Transition Analysis] Fate analysis failed: {fate_exc}")

    result['fate_result'] = fate_result
    result['df_condition_cmp'] = df_condition_comp
    result['df_paths'] = df_lineage_paths
    result['adata'] = a

    return result
