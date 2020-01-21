import os
import pandas as pd
import numpy as np

def volcano(file_path, inhibitor, p_val_threshold, fold_threshold): #!! change this, put in dif file
    df = pd.DataFrame(pd.read_csv(file_path, sep='\t'))
    inhibitor = inhibitor + '_'
    col_names = df.columns
    for name in col_names:
        if inhibitor in name:
            df.rename(columns = {name: name.replace(inhibitor,'')}, inplace = True)

    df.dropna(axis= 1, how='all', inplace = True)
    list_nulls = list(df[df.fold_change.isnull()].Substrate)
    list_control_zero = list(df[df.control_mean == 0].Substrate)
    df = df[(df.fold_change.notnull()) & (df.control_mean != 0) & (df.fold_change != 0)]
    df['log_foldchange'] = np.log2(df['fold_change'])
    df['log_pval'] = -np.log10(df['p-value'])

    df =df[~(np.isinf(np.abs(df.log_foldchange)))& ~(np.isinf(np.abs(df.log_pval))) & ~np.isnan(df.log_foldchange) & ~np.isnan(df.log_pval)]

    pval_threshold = int(p_val_threshold)
    fold_threshold = int(fold_threshold)
    x_min = df['log_foldchange'].min()
    x_max = df['log_foldchange'].max()
    y_max = df['log_pval'].max()
    plot_layout = {'x_min': x_min, 'x_max': x_max, 'y_max':y_max}

    df_negfold_pval = df[(df.log_pval > pval_threshold)  & (df.log_foldchange < -fold_threshold) ]
    df_posfold_pval = df[(df.log_pval > pval_threshold)  & (df.log_foldchange > fold_threshold) ]
    df_above_threshold = df[((df.log_pval <= pval_threshold)) | ((np.abs(df.log_foldchange) < fold_threshold)  )]

    out_dict = {}
    negfold_pval = {'substrate': list(df_negfold_pval.Substrate), 'foldchange': list(df_negfold_pval.log_foldchange), 'pval': list(df_negfold_pval.log_pval)}
    posfold_pval = {'substrate': list(df_posfold_pval.Substrate), 'foldchange': list(df_posfold_pval.log_foldchange), 'pval': list(df_posfold_pval.log_pval)}
    above_pval = {'substrate': list(df_above_threshold.Substrate), 'foldchange': list(df_above_threshold.log_foldchange), 'pval': list(df_above_threshold.log_pval)}
    out_dict['negfold_pval'] = negfold_pval
    out_dict['posfold_pval'] = posfold_pval
    out_dict['above_pval'] = above_pval
    out_dict['plot_layout'] = plot_layout
    os.remove(file_path)
    return out_dict
