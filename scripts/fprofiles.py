#!/usr/bin/env python3
"""
Script to perform Non-negative Matrix Factorization (NMF) on motif frequency data from a TSV file,
followed by deconvolution analysis using Non-Negative Least Squares (NNLS).
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.optimize import nnls
import logging
import argparse

# --- Load Arguments --- #
def load_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data_file", type=str, required=True, help="Path to the TSV file containing motif frequencies")
    parser.add_argument("--output_fprof_per_lib", type=str, required=True, help="Path to output file for F-profiles per library")
    parser.add_argument("--output_motif_per_fprof", type=str, required=True, help="Path to output file for normalized motif F-profiles")
    args = parser.parse_args()
    return args

# --- Helper Functions --- #
def load_data(filepath):
    """Load motif data from a TSV file."""
    return pd.read_csv(filepath, delimiter='\t', dtype=str)

def perform_nmf(data_matrix, n_components=6):
    """Perform NMF on the provided data matrix."""
    model = NMF(n_components=n_components, init='random', random_state=0, max_iter=100000)
    return model.fit_transform(data_matrix), model.components_

def perform_nnls(F, sample_motif_frequencies):
    """Perform NNLS to decompose sample motif frequencies into F-profile contributions."""
    P, _ = nnls(F.T, sample_motif_frequencies)
    return 100 * P / np.sum(P)

# --- Main Function --- #
def main():
    logging.basicConfig(level=logging.INFO)
    args = load_arguments()

    logging.info("Loading motif frequency data...")
    motifs_df = load_data(args.data_file)
    M = motifs_df.iloc[:, 1:].astype(float).T

    logging.info(f"Data loaded with shape {M.shape}. Performing NMF...")
    W, F = perform_nmf(M, n_components=6)  # Ensure n_components is defined or pass it as a parameter

    nnls_contributions_df = pd.DataFrame()
    for lib in motifs_df.columns[1:]:
        logging.info(f"Processing sample {lib}...")
        sample_motif_frequencies = motifs_df[lib].values.astype(float)
        P_normalized = perform_nnls(F, sample_motif_frequencies)
        nnls_contributions_df[lib] = P_normalized

    nnls_contributions_df.index = [f'fprof{i+1}' for i in range(F.shape[0])]
    nnls_contributions_df.to_csv(args.output_fprof_per_lib, sep='\t')
    logging.info(f"NNLS contributions saved to {args.output_fprof_per_lib}.")

    F_normalized = F / F.sum(axis=1, keepdims=True)
    F_df_normalized = pd.DataFrame(F_normalized, columns=motifs_df.iloc[:, 0], index=[f'fprof{i+1}' for i in range(F.shape[0])])
    F_df_normalized.to_csv(args.output_motif_per_fprof, sep='\t', index=True)
    logging.info(f"Normalized motif profiles saved to {args.output_motif_per_fprof}.")

# --- Main Guard --- #
if __name__ == "__main__":
    main()
