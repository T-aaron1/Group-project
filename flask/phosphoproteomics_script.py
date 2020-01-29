import os
import pandas as pd
import numpy as np
from scipy.stats import norm



def change_column_names(file_path, inhibitor): #!! change this, put in dif file
    df = pd.DataFrame(pd.read_csv(file_path, sep='\t'))
    os.remove(file_path)
    inhibitor = inhibitor.upper().rstrip() + '_'
    col_names = df.columns
    for name in col_names:
        if inhibitor in name:
            df.rename(columns = {name: name.replace(inhibitor,'')}, inplace = True)
    df.dropna(axis= 1, how='all', inplace = True)
    return df

#df instead of file_path...
def volcano(df, p_val_threshold, fold_threshold): #!! change this, put in dif file
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
    return out_dict

####

def extract_above_threshold(df, volcano_results):
    substrate = volcano_results['negfold_pval']['substrate'] + volcano_results['posfold_pval']['substrate']
    print(len(substrate))
    df_substrate = pd.DataFrame.from_dict({'Substrate':substrate})
    df2 = df_substrate.join(df.set_index(['Substrate']), on=['Substrate'])
    return df2

#Substrate  control_mean          mean  fold_change   p-value    ctrlCV   treatCV

#####

def KSEA(df, kinase_substrate):
    df.dropna(axis=1, how='all', inplace=True)

    df['subst_position'] = ''
    df['subst_position'] = df.Substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    df['subst_name'] = df.Substrate.str.split("(", n=1, expand=True)[0]

    df = df.dropna(axis=0)

    df = df.join(kinase_substrate[['kinase', 'sub_gene', 'sub_mod_rsd']].set_index(['sub_gene', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'])

    df = df.join(kinase_substrate[['kinase', 'sub_acc_id', 'sub_mod_rsd']].set_index(['sub_acc_id', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'], rsuffix='_sub_acc')

    df['kinase_final'] = df[['kinase', 'kinase_sub_acc']].fillna('').sum(axis=1)
    df = df[df.fold_change != 'inf']
    df['Log2Substrate_fold_change'] = np.log2(df.fold_change)
    mean_log2_FC = df["Log2Substrate_fold_change"].mean()
    kinase_count = df['kinase_final'].value_counts()
    kinase_count_df = pd.DataFrame(kinase_count)
    kinase_count_df['sqrt_kinase'] = np.sqrt(kinase_count_df)
    df.groupby('kinase_final')['Log2Substrate_fold_change'].mean()
    kinase_count_df['Standard_deviation'] = np.std(df.Log2Substrate_fold_change)
    kinase_count_df['Log2Substrate_fold_change'] = df.groupby('kinase_final')['Log2Substrate_fold_change'].mean()
    kinase_count_df['mean_log2(FC)'] = df['mean_log2(FC)'] = df["Log2Substrate_fold_change"].mean()
    kinase_count_df['KSEA'] = (kinase_count_df['Log2Substrate_fold_change'] - mean_log2_FC *
                               kinase_count_df['sqrt_kinase']) / kinase_count_df['Standard_deviation']


    kinase_count_df['P_value'] = norm.sf(abs(kinase_count_df['KSEA'])) * 2
    non_identified = kinase_count_df.loc['', 'kinase_final']

    kinase_count_df.drop([''], axis=0, inplace=True)
    kinase_count_df = kinase_count_df.sort_values(['KSEA'], ascending=True)
    kinase_count_df['Kinase'] = kinase_count_df.index
    kinase_count_df.reset_index(drop=True, inplace=True)
    sig_kinase_count = kinase_count_df[kinase_count_df['P_value']<0.05]
    colors = ["rgb(255,242,0)"]*sig_kinase_count[sig_kinase_count['KSEA']<0].shape[0] +  ["rgb(0,34,255)"]*sig_kinase_count[sig_kinase_count['KSEA']>=0].shape[0]
    kinase_count_dict={'kinase': list(sig_kinase_count.loc[:,'Kinase']),
                       'ksea' : list(sig_kinase_count.loc[:,'KSEA']),
                       'p_value': list(sig_kinase_count.loc[:,'P_value']),
                       'colors':colors}

    output = {'score': kinase_count_dict ,'non_identified': non_identified}

    return output
