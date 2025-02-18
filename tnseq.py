import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import sys

def read_wig(file, window_size):
    df = pd.read_csv(file, sep='\t', header=None, names=['Position', 'Count', 'Gene'])

    # Split 'Gene' column by comma and expand into separate rows
    df = df.assign(Gene=df['Gene'].str.split(', ')).explode('Gene')

    # Remove positions with counts at a frequency greater than 10^-3
    total_counts = df['Count'].sum()
    df = df[df['Count'] / total_counts <= 10**-3]

    # Normalize by scaling all counts to 10 million reads
    df['Count'] = df['Count'] * (10**7 / df['Count'].sum())

    # Sum counts in the determined window size for gene positions only
    df_gene = df.dropna().copy()  # create a copy of the DataFrame to avoid SettingWithCopyWarning
    df_gene.loc[:, 'Window'] = df_gene['Position'] // window_size
    df_gene = df_gene.groupby(['Gene', 'Window'])['Count'].sum().reset_index()

    return df_gene.groupby('Gene')['Count'].apply(list)

def calculate_window_size(file):
    df = pd.read_csv(file, sep='\t', header=None, names=['Position', 'Count', 'Gene'])
    df['Gene'] = df['Gene'].fillna(method='ffill')

    # Calculate the frequency of non-zero positions
    non_zero_positions = df[df['Count'] != 0].shape[0]
    total_positions = df.shape[0]
    frequency = non_zero_positions / total_positions

    # Determine the window size
    window_size = int(1 / frequency)

    return window_size, non_zero_positions

def log2_fold_change(control, experiment):
    return np.log2(np.mean(experiment) / np.mean(control))

def mann_whitney_u_test(control, experiment):
    _, p_value = mannwhitneyu(control, experiment)
    return p_value

def main(control_file, experiment_file):
    control_window_size, control_non_zero_positions = calculate_window_size(control_file)
    experiment_window_size, experiment_non_zero_positions = calculate_window_size(experiment_file)
    window_size = max(control_window_size, experiment_window_size)

    control = read_wig(control_file, window_size)
    experiment = read_wig(experiment_file, window_size)
    all_genes = control.index.union(experiment.index)
    control = control.reindex(all_genes, fill_value=[0])
    experiment = experiment.reindex(all_genes, fill_value=[0])
    genes = pd.DataFrame(index=all_genes)
    genes['Log2FoldChange'] = [log2_fold_change(control[gene], experiment[gene]) for gene in genes.index]
    p_values = [mann_whitney_u_test(control[gene], experiment[gene]) for gene in genes.index]
    genes['P-value'] = p_values

    # Adjust p-values for multiple testing using Benjamini-Hochberg procedure
    _, p_values_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
    genes['Adjusted P-value'] = p_values_adjusted
    genes['MinusLog10(Adjusted P-value)'] = -np.log10(p_values_adjusted)

    # Replace -inf in Log2FoldChange with the next lowest whole integer less than the smallest non-infinite value
    min_value = genes[genes['Log2FoldChange'] != -np.inf]['Log2FoldChange'].min()
    genes['Log2FoldChange'] = genes['Log2FoldChange'].replace(-np.inf, np.floor(min_value) - 0)

    # Create output file name
    control_name = control_file.split('_')[0]
    experiment_name = experiment_file.split('_')[0]
    output_file = f"{control_name}_{experiment_name}_mwu_analysis.csv"

    genes.to_csv(output_file)

    print(f"Total non-zero positions in control file: {control_non_zero_positions}")
    print(f"Total non-zero positions in experiment file: {experiment_non_zero_positions}")
    print(f"Calculated window size: {window_size}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])