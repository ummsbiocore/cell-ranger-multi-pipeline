#!/usr/bin/env python3

import argparse
import sys
import time
import numpy as np
import pandas as pd
from pathlib import PurePath

# Try to import regdiffusion
try:
    import regdiffusion as rd
except ImportError:
    print("Error: regdiffusion not found. Please install it in your environment.", file=sys.stderr)
    sys.exit(1)

# Try to import pyscenic utils
try:
    from pyscenic.cli.utils import load_exp_matrix, suffixes_to_separator
    from arboreto.utils import load_tf_names
except ImportError:
    print("Error: pyscenic not found. Please install it in your environment.", file=sys.stderr)
    sys.exit(1)

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description="Run RegDiffusion (GPU accelerated) for GRN inference, as an alternative to Arboreto."
    )

    parser.add_argument(
        "expression_mtx_fname",
        type=str,
        help="The name of the file that contains the expression matrix for the single cell experiment."
        " Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).",
    )
    parser.add_argument(
        "tfs_fname",
        type=str,
        help="The name of the file that contains the list of transcription factors (TXT; one TF per line).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=sys.stdout,
        help="Output file/stream, i.e. a table of TF-target genes (TSV).",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=4,
        help="The number of workers to use for edge extraction. (default: 4).",
    )
    parser.add_argument(
        "--cell_id_attribute",
        type=str,
        default="CellID",
        help="The name of the column attribute that specifies the identifiers of the cells in the loom file.",
    )
    parser.add_argument(
        "--gene_attribute",
        type=str,
        default="Gene",
        help="The name of the row attribute that specifies the gene symbols in the loom file.",
    )
    parser.add_argument(
        "--sparse",
        action="store_const",
        const=True,
        default=False,
        help="If set, load the expression data as a sparse (CSC) matrix.",
    )
    parser.add_argument(
        "-t",
        "--transpose",
        action="store_const",
        const="yes",
        help="Transpose the expression matrix (rows=genes x columns=cells).",
    )
    parser.add_argument(
        "--no_log_transform",
        action="store_true",
        help="If set, skip the log(x+1) transformation. Use this if your data is already log-transformed.",
    )
    parser.add_argument(
        "--top_gene_percentile",
        type=float,
        default=50.0,
        help="Percentile threshold for keeping edges based on weights (default: 50.0).",
    )
    parser.add_argument(
        "--num_hvg",
        type=int,
        default=0,
        help="Number of Highly Variable Genes (HVGs) to select. Set to 0 to use all genes. (default: 0)",
    )
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Run on CPU instead of GPU.",
    )

    return parser

def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    start_time = time.time()
    
    # Load expression matrix
    # load_exp_matrix returns (matrix, gene_names) if sparse, else DataFrame
    ex_matrix_obj = load_exp_matrix(
        args.expression_mtx_fname,
        (args.transpose == "yes"),
        args.sparse,
        args.cell_id_attribute,
        args.gene_attribute,
    )

    if args.sparse:
        gene_names = ex_matrix_obj[1]
        ex_matrix = ex_matrix_obj[0]
        # Convert to dense for RegDiffusion
        # Warning: This might consume a lot of memory for very large datasets
        print("Converting sparse matrix to dense for RegDiffusion...", file=sys.stdout)
        ex_matrix = ex_matrix.toarray()
    else:
        gene_names = ex_matrix_obj.columns
        ex_matrix = ex_matrix_obj.values

    end_time = time.time()
    print(
        f"Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time:.2f} seconds...",
        file=sys.stdout,
    )

    # Load TFs
    tf_names = load_tf_names(args.tfs_fname)
    print(f"Loaded {len(tf_names)} TFs...", file=sys.stdout)

    # Preprocessing
    if not args.no_log_transform:
        print("Applying log(x+1) transformation...", file=sys.stdout)
        ex_matrix = np.log(ex_matrix + 1.0)
    else:
        print("Skipping log transformation as requested...", file=sys.stdout)

    # HVG Selection
    if args.num_hvg > 0 and args.num_hvg < len(gene_names):
        print(f"Selecting top {args.num_hvg} Highly Variable Genes (HVGs)...", file=sys.stdout)
        try:
            import scanpy as sc
            import anndata as ad
            
            # Create AnnData object
            adata = ad.AnnData(X=ex_matrix, var=pd.DataFrame(index=gene_names))
            
            # Calculate HVGs
            # We assume data is already log-transformed (or user said no_log_transform)
            # flavor='seurat' expects log-transformed data
            sc.pp.highly_variable_genes(adata, n_top_genes=args.num_hvg, subset=True)
            
            # Update matrix and gene names
            ex_matrix = adata.X
            gene_names = adata.var_names
            print(f"Data shape after HVG filtering: {ex_matrix.shape}", file=sys.stdout)
            
        except ImportError:
            print("Warning: scanpy or anndata not found. Skipping HVG selection.", file=sys.stderr)
        except Exception as e:
            print(f"Warning: HVG selection failed: {e}. Using all genes.", file=sys.stderr)

    # Determine device
    device = 'cpu' if args.cpu else 'cuda'
    print(f"Using device: {device}", file=sys.stdout)

    # Run RegDiffusion
    print("Initializing RegDiffusionTrainer...", file=sys.stdout)
    rd_trainer = rd.RegDiffusionTrainer(ex_matrix, device=device)
    
    print(f"Starting training on {device}...", file=sys.stdout)
    train_start = time.time()
    rd_trainer.train()
    print(f"Training finished in {time.time() - train_start:.2f} seconds.", file=sys.stdout)

    # Get GRN
    print("Extracting GRN...", file=sys.stdout)
    # Using top_gene_percentile as suggested in tutorial/docs to prune weak edges initially if needed
    # But to mimic arboreto we might want more edges. 
    # However, RegDiffusion tutorial suggests filtering.
    grn = rd_trainer.get_grn(gene_names, top_gene_percentile=args.top_gene_percentile)

    # Extract edge list
    # k=-1 to extract all edges (from the kept percentile), then we filter by TFs
    print(f"Extracting edge list with {args.num_workers} workers...", file=sys.stdout)
    edgelist = grn.extract_edgelist(k=-1, workers=args.num_workers)
    
    # Rename columns to match pySCENIC expectations if needed, or just standard TF/target/importance
    # RegDiffusion extract_edgelist returns columns: ['regulator', 'target', 'weight'] usually?
    # Tutorial says: edgelist.columns = ['TF', 'target', 'importance']
    edgelist.columns = ['TF', 'target', 'importance']

    # Filter for TFs in our list
    print("Filtering for provided Transcription Factors...", file=sys.stdout)
    edgelist = edgelist[edgelist['TF'].isin(tf_names)]

    # Sort by importance
    edgelist = edgelist.sort_values(by="importance", ascending=False)

    end_time = time.time()
    print(f"Done in {end_time - start_time:.2f} seconds.", file=sys.stdout)

    # Save output
    if args.output == sys.stdout:
        # If writing to stdout, just print csv
        edgelist.to_csv(sys.stdout, index=False, sep='\t')
    else:
        extension = PurePath(args.output).suffixes
        sep = suffixes_to_separator(extension)
        print(f"Saving results to {args.output}...", file=sys.stdout)
        edgelist.to_csv(args.output, index=False, sep=sep)

if __name__ == "__main__":
    main()
