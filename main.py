import os
import math 
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.stats as stats

from scipy.cluster.hierarchy import linkage
from plotly.figure_factory import create_dendrogram
from plotly.subplots import make_subplots
from statsmodels.stats.multitest import multipletests


def get_df_dv_per_region(name, ds, regions, adata, debug, p_lim):
    all_results = []

    for region in regions:
        print(f"Processing {region} for Region # {name}...")
        if debug:
            print(f"Processing region: {region}")
        # Get differential velocity results for the current region
        region_results = get_diff_velo_genes(name, ds, region, adata, debug, p_lim)
        
        if region_results is not None and not region_results.empty:
            # Add the region name to the results
            region_results['region'] = region
            all_results.append(region_results)

        print(f"{region} for {name}: done.")

    # Combine all region results into a single DataFrame
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
    else:
        combined_results = pd.DataFrame()  # Return an empty DataFrame if no results
    
    return combined_results


def get_diff_velo_genes(name, ds, region, adata, debug = False, p_lim=1e-9, plot=True, save_plot_path='difference_distribution.png'):
    """
    Identify differentially veloce genes
    
    Parameters:
    -----------
    ds: int
    region : str
        Region of interest for analysis
    adata : AnnData
        Annotated data matrix
    p_lim : float, default 1e-9
        P-value threshold for significance
    
    Returns:
    --------
    pd.DataFrame
        DataFrame of significant differentially velocized genes
    """
    # Extract velocity data
    velocity_df = pd.DataFrame(
        adata.layers['velocity'],
        index=adata.obs.index,
        columns=adata.var.index
    )
    velocity_df_normalized = velocity_df # (velocity_df - velocity_df.mean(axis=0)) / velocity_df.std(axis=0)    
    #print(velocity_df_normalized)
       
    # Add condition and region information
    velocity_df_normalized['Condition'] = adata.obs['Condition']
    velocity_df_normalized['isCPU'] = adata.obs['Regions'] == region

    # Filter for specific region
    velocity_df_n_CPU = velocity_df_normalized[velocity_df_normalized['isCPU']]
    velocity_df_n_CPU = velocity_df_n_CPU.drop('isCPU', axis=1)
    
    # Predefined conditions (consider making these parameters)
    condition_1 = ['CF', 'CM', 'CF','CM']
    condition_2 = ['EcigF', 'EcigM', 'EcigF', 'EcigM']
    
    # Separate conditions
    velocity_group_1 = velocity_df_n_CPU[velocity_df_n_CPU['Condition'] == condition_1[ds]].drop('Condition', axis=1)
    velocity_group_2 = velocity_df_n_CPU[velocity_df_n_CPU['Condition'] == condition_2[ds]].drop('Condition', axis=1)

    results = []
    if debug:
        print(velocity_df_n_CPU)
    try:
        for gene in velocity_group_1.columns:
            if gene not in velocity_group_2.columns:
                continue
                
            # Drop NaN values for the current gene in both conditions
            group_1_valid = velocity_group_1[gene].dropna()
            group_2_valid = velocity_group_2[gene].dropna()

            # Proceed only if both groups have sufficient data
            if len(group_1_valid) > 1 and len(group_2_valid) > 1:
                # Use Welch's t-test (does not assume equal variances)
                stat, pval = stats.ttest_ind(
                    group_1_valid, 
                    group_2_valid, 
                    nan_policy='omit', 
                    equal_var=False
                )
                
                # Additional metrics for vector normalization
                diff = math.exp(np.mean(group_1_valid)) / (math.exp(np.mean(group_1_valid)) + math.exp(np.mean(group_2_valid)))
                #print(diff)
                #return None
                results.append({
                    'gene': gene, 
                    't_statistic': stat, 
                    'p_value': pval,
                    'size_control': len(group_1_valid),
                    'size_treated': len(group_2_valid),
                    'difference': diff,
                })
            else:
                print(f"Empty group for: {gene}")
                # Append a placeholder for insufficient data
                results.append({
                    'gene': gene, 
                    't_statistic': None, 
                    'p_value': None,
                    'size_control': None,
                    'size_treated': None,                    
                    'difference': None
                })
    except Exception as e:
        print(f"Error during gene analysis: {e}")
        return None
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    if debug:
        print(results_df)
    # Remove rows with None or NaN p-values
    results_df = results_df.dropna(subset=['p_value'])

    try:
        # Adjust for multiple testing (Benjamini-Hochberg)
        results_df['adjusted_p_value'] = multipletests(
            results_df['p_value'], 
            method='fdr_bh'
        )[1]
        
        # Filter significant genes
        significant_genes = results_df[results_df['adjusted_p_value'] <= p_lim]
        save_plot_path = name + "_" + region.replace("/","-") + "_" + save_plot_path
        if plot:
            plot_difference_distribution(
                results_df['difference'],
                save_path=save_plot_path,
                title='Difference Distribution of Genes',
                xlabel='Difference',
                ylabel='Density'
            )
        return significant_genes
    
    except Exception as e:
        print(f"Error during multiple testing correction: {e}")


def generate_interactive_clustermap_with_dendrogram(
    ds,
    df,
    sub_folder_name,
    value_col='adjusted_p_value',
    gene_col='gene',
    region_col='region',
    output_dir='output',
    cmap='Viridis',
    width=1200,
    height=800
):
    """
    Generates an interactive cluster heatmap with visible dendrograms for rows (genes)
    and columns (regions). Debugging outputs are added to check for issues in the pivot
    table and data preparation steps.

    Steps:
    1) Pivot df into gene x region matrix
    2) Replace small / NaN values
    3) -log10 transform
    4) Create row dendrogram -> get ordering of genes
    5) Create column dendrogram -> get ordering of regions
    6) Reorder matrix
    7) Build a subplot layout:
       [ top dendrogram  ][  top dendrogram ]
       [ row dendrogram  ][    heatmap     ]
    8) Make each plot interactive, so you can hover over dendrogram lines and heatmap cells
    9) Save an HTML file; user can open it in any browser

    :param ds: str, dataset name
    :param df: pd.DataFrame, must have (gene_col, region_col, value_col)
    :param sub_folder_name: str, subdirectory inside output_dir for saving results
    :param value_col: str, name of values column in df (default 'adjusted_p_value')
    :param gene_col: str, name of gene column (default 'gene')
    :param region_col: str, name of region column (default 'region')
    :param output_dir: str, parent directory for saving
    :param cmap: str, Plotly color scale
    :param width: int, width of the overall figure in pixels
    :param height: int, height of the overall figure in pixels
    """

    # 1. Create output folders
    os.makedirs(output_dir, exist_ok=True)
    sub_folder_path = os.path.join(output_dir, sub_folder_name)
    os.makedirs(sub_folder_path, exist_ok=True)

    # 2. Debug: Check unique regions and their counts
    unique_regions = df[region_col].unique()
    print(f"Found {len(unique_regions)} unique regions: {unique_regions}")
    print(f"Region counts:\n{df[region_col].value_counts()}")

    # 3. Pivot the DataFrame
    try:
        matrix = df.pivot_table(
            values=value_col,
            index=gene_col,
            columns=region_col,
            aggfunc='min'
        )
    except Exception as e:
        print("Error during pivot_table:", e)
        return

    print("After pivot:")
    print("  shape =", matrix.shape)
    print("  columns =", matrix.columns.tolist())

    # 4. Save the input df to CSV for reference
    csv_path = os.path.join(sub_folder_path, f'{ds}_heatmap_data.csv')
    df.to_csv(csv_path, index=False)
    print(f"Input data saved to {csv_path}")

    # 5. Replace extremely small / zero values, fill NaNs
    min_val = 1e-100
    matrix[matrix < min_val] = min_val
    matrix.fillna(1.0, inplace=True)

    print("Matrix after replacing small/NaN values:")
    print(matrix.describe())

    # 6. -log10 transform
    matrix = np.log10(matrix)

    print("Matrix after -log10 transformation:")
    print(matrix.describe())

    # 7. Create row dendrogram
    row_dendro = create_dendrogram(
        matrix.values,
        orientation='right',
        linkagefun=lambda x: linkage(x, 'ward', metric='euclidean')
    )
    row_leaves = row_dendro['layout']['yaxis']['ticktext']
    row_leaves_idx = list(map(int, row_leaves))

    # 8. Create column dendrogram
    col_dendro = create_dendrogram(
        matrix.values.T,
        orientation='bottom',
        linkagefun=lambda x: linkage(x, 'ward', metric='euclidean')
    )
    col_leaves = col_dendro['layout']['xaxis']['ticktext']
    col_leaves_idx = list(map(int, col_leaves))

    # 9. Reorder the matrix
    matrix_reordered = matrix.iloc[row_leaves_idx, col_leaves_idx]
    reordered_row_labels = matrix.index[row_leaves_idx]
    reordered_col_labels = matrix.columns[col_leaves_idx]

    print("Reordered matrix shape:", matrix_reordered.shape)
    print("Reordered column labels:", reordered_col_labels)
    print("Reordered row labels:", reordered_row_labels[:10])  # Only print first 10 for brevity

    # 10. Build composite figure with subplots
    fig = make_subplots(
        rows=2, cols=2,
        row_heights=[0.2, 0.8],
        column_widths=[0.2, 0.8],
        horizontal_spacing=0.02,
        vertical_spacing=0.02,
        specs=[
            [{"type": "xy", "rowspan": 1, "colspan": 1}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}]
        ]
    )

    # Add column dendrogram (top, row=1 col=2)
    for trace in col_dendro['data']:
        fig.add_trace(trace, row=1, col=2)
    fig.update_xaxes(autorange='reversed', row=1, col=2)
    fig.update_xaxes(showticklabels=False, row=1, col=2)

    # Add row dendrogram (left, row=2 col=1)
    for trace in row_dendro['data']:
        fig.add_trace(trace, row=2, col=1)
    fig.update_yaxes(showticklabels=False, row=2, col=1)

    # Add heatmap (row=2 col=2)
    heatmap = go.Heatmap(
        z=matrix_reordered.values,
        x=reordered_col_labels,
        y=reordered_row_labels,
        colorscale=cmap,
        hovertemplate="Gene: %{y}<br>Region: %{x}<br>Value: %{z}<extra></extra>"
    )
    fig.add_trace(heatmap, row=2, col=2)

    fig.update_xaxes(side="top", row=2, col=2, tickangle=45)
    fig.update_yaxes(autorange='reversed', row=2, col=2)

    fig.update_layout(
        width=width,
        height=height,
        title_text=f"Interactive Clustered Heatmap & Dendrogram: {ds.upper()}"
    )

    # 11. Save HTML
    html_path = os.path.join(sub_folder_path, f'{ds}_interactive_clustermap.html')
    fig.write_html(html_path)
    print(f"Interactive dendrogram heatmap saved to {html_path}")

    # 12. Show
    fig.show()

def get_adata(idx):
    FILE_ADATA = [
        "<PATH>.file_1.h5ad",
        "<PATH>.file_2.h5ad",
        "<PATH>.file_3.h5ad",
        "<PATH>.file_4.h5ad",
        "<PATH>.file_5.h5ad",
        "..."
    ]
    DATASETS = ["DS1","DS2","DS3","DS4"]
    
    intg_adata_gv = sc.read_h5ad(FILE_ADATA[idx])
    return DATASETS[idx], intg_adata_gv


def run_analysis(plim_lowest, plim_name):
    for i in range(0,4):
        try:
            ds, adata = get_adata(i)
            tmp_r_names = [str(region) for region in adata.obs['Regions'].unique().tolist()]
            r_names = []
            for region in tmp_r_names:
                if region not in ["Layer1", "LV", "Sparse"]:
                    r_names.append(region)
                    
            df = get_df_dv_per_region(ds, i, r_names, adata, False, plim_lowest)
            generate_interactive_clustermap_with_dendrogram(ds, df, plim_name)
        except Exception as e:
            print(e)


run_analysis(1e-18, "final")
