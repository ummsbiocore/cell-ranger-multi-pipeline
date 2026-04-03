#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import argparse
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os
import sys
import pickle
from analysis_utils import perform_transition_analysis, run_velocity_pipeline

def main():
    parser = argparse.ArgumentParser(description="Pre-compute scVelo transition analysis results.")
    parser.add_argument("input_file", help="Path to input .h5ad file")
    parser.add_argument("--output_file", help="Path to output .h5ad file (default: overwrite input)", default=None)
    parser.add_argument("--group_col", required=True, help="Column for clustering/grouping")
    parser.add_argument("--condition_col", help="Column for experimental conditions")
    parser.add_argument("--k_steps", default="1,2,3", help="Comma-separated list of k steps (default: 1,2,3)")
    parser.add_argument("--n_macrostates", type=int, default=10, help="Number of macrostates (default: 10)")
    parser.add_argument("--target_clusters", help="Comma-separated list of target clusters to filter (optional)")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of threads (default: 1)")

    args = parser.parse_args()

    input_path = args.input_file
    output_path = args.output_file if args.output_file else input_path
    
    # Set random seed for reproducibility
    np.random.seed(0)
    scv.settings.seed = 0
    
    print(f"Loading data from {input_path}...")
    adata = sc.read_h5ad(input_path)
    
    # Parse list arguments
    k_list = [int(k.strip()) for k in args.k_steps.split(',') if k.strip()]
    filter_set = set([c.strip() for c in args.target_clusters.split(',') if c.strip()]) if args.target_clusters else None
    
    # Determine threads
    n_jobs = args.n_jobs
    
    print(f"Using {n_jobs} threads.")

    # 1. Analyze All Cells
    print("Running transition analysis for All Cells...")
    global_result = perform_transition_analysis(
        adata,
        args.group_col,
        k_list,
        filter_set,
        args.n_macrostates,
        condition_col=args.condition_col,
        compute_condition_cmp=bool(args.condition_col),
        n_jobs=n_jobs
    )
    
    # The global_result['adata'] is the modified adata with velocity graph etc.
    # We should update our main adata with this
    # The global_result['adata'] is the modified adata with velocity graph etc.
    # We will use this for the final output, but keep 'adata' clean for condition analysis
    adata_out = global_result['adata']
    
    # Save lineage names and colors to uns to persist them through h5ad saving
    # Save lineage names and colors to uns to persist them through h5ad saving
    if 'lineages_fwd' in adata_out.obsm:
        if hasattr(adata_out.obsm['lineages_fwd'], 'names'):
            adata_out.uns['lineages_fwd_names'] = list(adata_out.obsm['lineages_fwd'].names)
        if hasattr(adata_out.obsm['lineages_fwd'], 'colors'):
            adata_out.uns['lineages_fwd_colors'] = list(adata_out.obsm['lineages_fwd'].colors)
            
    # Compute PAGA layout and save it
    # analysis_utils runs scv.tl.paga (graph), but we need the embedding positions
    # Compute PAGA layout and save it
    # analysis_utils runs scv.tl.paga (graph), but we need the embedding positions
    if 'paga' in adata_out.uns:
        print("Computing PAGA layout for consistent visualization...")
        # scv.pl.paga computes 'pos' in uns['paga']
        # We use a dummy plot call or just compute layout if possible.
        # scv.pl.paga calls scv.tools.paga if needed, but we already have it.
        # It calls scv.utils.get_paga_embedding -> sc.tl.draw_graph or similar.
        # Let's just call scv.pl.paga with show=False to generate the positions.
        try:
            scv.pl.paga(adata_out, show=False)
            print("PAGA layout computed and saved.")
        except Exception as e:
            print(f"Warning: Failed to compute PAGA layout: {e}")

    # Prepare results dictionary for caching
    # We remove 'adata' from the dictionary to avoid recursion/bloat when saving to uns
    global_result_clean = {k: v for k, v in global_result.items() if k != 'adata'}
    
    condition_results = {}
    
    # 2. Analyze Conditions
    if args.condition_col and args.condition_col in adata.obs.columns:
        conditions = adata.obs[args.condition_col].unique()
        print(f"Found conditions: {conditions}")
        
        for cond in conditions:
            print(f"Analyzing condition: {cond}")
            subset = adata[adata.obs[args.condition_col] == cond].copy()
            print(f"Subset size for {cond}: {subset.n_obs} cells")
            
            try:
                # Determine threads for subset
                subset_n_jobs = n_jobs
                
                # Recompute velocity graph and neighbors for the subset
                # This is crucial to avoid dimension mismatches and corrupted graphs
                print(f"Recomputing velocity graph for condition: {cond}")
                subset = run_velocity_pipeline(subset, n_jobs=subset_n_jobs)
                
                cond_res = perform_transition_analysis(
                    subset,
                    args.group_col,
                    k_list,
                    filter_set,
                    args.n_macrostates,
                    condition_col=None,
                    compute_condition_cmp=False,
                    n_jobs=subset_n_jobs
                )
                
                # Store heavy arrays in main adata with suffix
                # We need to map subset indices to main adata indices?
                # Actually, simply storing them in uns/obsm with suffix is tricky because shapes differ.
                # But wait, if we use the subset adata in the app, we usually recompute.
                # If we want to avoid recomputation, we need to save the results.
                # The most important results are:
                # - lineages_fwd (obsm)
                # - velocity_graph (uns - sparse matrix)
                # - velocity_pseudotime (obs)
                # - paga (uns) - NEW: Save PAGA results
                
                # Strategy: Store these in a dictionary in adata.uns['precomputed_condition_data'][cond]
                # We can store the sparse matrix and arrays there.
                # When loading in app, we inject them into the subset adata.
                
                cond_data = {}
                subset_adata = cond_res['adata']
                
                if 'lineages_fwd' in subset_adata.obsm:
                    cond_data['lineages_fwd'] = subset_adata.obsm['lineages_fwd']
                    # Also store names/colors if present
                    if hasattr(subset_adata.obsm['lineages_fwd'], 'names'):
                        cond_data['lineages_fwd_names'] = list(subset_adata.obsm['lineages_fwd'].names)
                    if hasattr(subset_adata.obsm['lineages_fwd'], 'colors'):
                        cond_data['lineages_fwd_colors'] = list(subset_adata.obsm['lineages_fwd'].colors)

                if 'velocity_graph' in subset_adata.uns:
                    cond_data['velocity_graph'] = subset_adata.uns['velocity_graph']
                    
                if 'velocity_graph_neg' in subset_adata.uns:
                    cond_data['velocity_graph_neg'] = subset_adata.uns['velocity_graph_neg']

                if 'velocity_pseudotime' in subset_adata.obs.columns:
                    # Store as Series with index to ensure alignment
                    cond_data['velocity_pseudotime'] = subset_adata.obs['velocity_pseudotime']
                    
                if 'paga' in subset_adata.uns:
                    cond_data['paga'] = subset_adata.uns['paga']

                # Clean result dict
                cond_res_clean = {k: v for k, v in cond_res.items() if k != 'adata'}
                cond_res_clean['precomputed_data'] = cond_data
                
                condition_results[str(cond)] = cond_res_clean
                
                # Explicitly clear memory
                del subset
                del cond_res
                del subset_adata
                import gc
                gc.collect()
                
            except Exception as e:
                print(f"Error analyzing condition {cond}: {e}")
                condition_results[str(cond)] = {'error': str(e)}

    # Construct final cache structure
    # Key format: f"{condition_col}_{group_col}_{k_list}_{n_macrostates}_{filter_vals}_{n_jobs}"
    # Note: n_jobs in key might be tricky if it varies. The app uses the current n_jobs to generate key.
    # If we precompute with n_jobs=X, and app runs with n_jobs=Y, key won't match.
    # However, the RESULTS shouldn't depend on n_jobs (it's just for speed).
    # So we should probably make the app ignore n_jobs in cache key, OR we force the app to use the precomputed key.
    # Better: Store in a fixed location `adata.uns['precomputed_results']` and have app check that first, ignoring the cache key check for n_jobs.
    
    # Actually, let's just use the parameters to build the key, but maybe store it in a way the app can find.
    # The app generates a key.
    # Let's store it in `adata.uns['precomputed_results']` as a dictionary where keys are the cache keys.
    
    # We need to match the app's key generation exactly.
    # App key: f"{condition_col}_{group_col}_{','.join(map(str, k_list))}_{n_macrostates}_{','.join(sorted(filter_vals)) if filter_vals else 'all'}_{n_jobs}"
    
    # We will assume the user will run the app with the same n_jobs or we should relax the key check in the app.
    # Let's relax the key check in the app to ignore n_jobs if a precomputed result is found.
    # For now, let's generate the key with the n_jobs we used.
    
    filter_str = ','.join(sorted(filter_set)) if filter_set else 'all'
    cache_key = f"{args.condition_col if args.condition_col else 'None'}_{args.group_col}_{','.join(map(str, k_list))}_{args.n_macrostates}_{filter_str}_{n_jobs}"
    
    full_result = dict(global_result_clean)
    full_result['condition_col'] = args.condition_col
    full_result['condition_results'] = condition_results
    
    if 'precomputed_results' not in adata_out.uns:
        adata_out.uns['precomputed_results'] = {}
    
    # We need to use a distinct key or structure.
    # Since uns keys must be strings, and values usually simple.
    # Storing a complex dict with DataFrames and sparse matrices in uns might be an issue for h5ad saving if not careful.
    # AnnData handles dicts in uns, but deep nesting of complex types can be problematic.
    # It's safer to pickle the result and store as a byte string, OR rely on AnnData's native support.
    # Let's try native support first, but if it fails, we might need to pickle.
    # Actually, pickling is safest for arbitrary python objects like our result dict.
    # But h5ad doesn't support bytes directly in uns easily without wrapping.
    # Let's try to store it as a special encoded string or just rely on the fact that we are in python.
    # Wait, if we save to h5ad, it must be hdf5 compatible.
    # DataFrames in uns are supported? Not really standard.
    # Usually DataFrames are in obs/var/obsm.
    # Our result has DataFrames (`cluster_matrix_by_k`, `df_condition_cmp`).
    
    # Alternative: Store the heavy data (matrices) in obsm/uns of adata, and the lightweight metadata in a JSON string.
    # But we have many matrices (one per k).
    
    # Let's use pickle and store as a void array or similar?
    # Or just use a separate file? User asked for "output h5ad file which includes results".
    
    # Let's try to structure it compatible with AnnData.
    # Convert DataFrames to dicts?
    # Or just use pickle.
    # To store pickle in h5ad:
    # import numpy as np
    # adata.uns['precomputed_blob'] = np.void(pickle.dumps(full_result))
    # This works for saving/loading in Python.
    
    print("Saving results to adata.uns['precomputed_results']...")
    # We'll use a simplified key for the precomputed result to make it easy to find
    # We'll store a list of precomputed configurations.
    
    # Actually, let's just pickle the whole result dict and store it.
    # It's the most robust way to preserve the exact structure including DataFrames.
    blob = np.frombuffer(pickle.dumps(full_result), dtype=np.uint8)
    
    # Store in a dict keyed by the parameters (excluding n_jobs to be flexible)
    # Key: f"{condition_col}_{group_col}_{k_list}_{n_macrostates}_{filter_str}"
    param_key = f"{args.condition_col if args.condition_col else 'None'}_{args.group_col}_{','.join(map(str, k_list))}_{args.n_macrostates}_{filter_str}"
    
    if 'precomputed_blobs' not in adata_out.uns:
        adata_out.uns['precomputed_blobs'] = {}
        
    adata_out.uns['precomputed_blobs'][param_key] = blob
    
    print(f"Saving modified AnnData to {output_path}...")
    adata_out.write_h5ad(output_path)
    print("Done.")
    os._exit(0)

if __name__ == "__main__":
    main()
