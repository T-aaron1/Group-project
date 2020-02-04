import os
import pandas as pd
import numpy as np
from scipy.stats import norm



def change_column_names(file_path, inhibitor): #!! change this, put in dif file
    df = pd.DataFrame(pd.read_csv(file_path, sep='\t'))
    os.remove(file_path)
    inhibitor = inhibitor.lower().rstrip() + '_'
    print(inhibitor)
    for name in df.columns:
        df.rename(columns = {name: name.lower()}, inplace = True)
    for name in df.columns:
        if inhibitor in name.lower():
            df.rename(columns = {name: name.lower().replace(inhibitor,'')}, inplace = True)
    df.dropna(axis= 1, how='all', inplace = True)
    return df

#df instead of file_path...
def volcano(df, p_val_threshold, fold_threshold): #!! change this, put in dif file
    list_nulls = list(df[df.fold_change.isnull()].substrate)
    list_control_zero = list(df[df.control_mean == 0].substrate)
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
    negfold_pval = {'substrate': list(df_negfold_pval.substrate), 'foldchange': list(df_negfold_pval.log_foldchange), 'pval': list(df_negfold_pval.log_pval)}
    posfold_pval = {'substrate': list(df_posfold_pval.substrate), 'foldchange': list(df_posfold_pval.log_foldchange), 'pval': list(df_posfold_pval.log_pval)}
    above_pval = {'substrate': list(df_above_threshold.substrate), 'foldchange': list(df_above_threshold.log_foldchange), 'pval': list(df_above_threshold.log_pval)}
    out_dict['negfold_pval'] = negfold_pval
    out_dict['posfold_pval'] = posfold_pval
    out_dict['above_pval'] = above_pval
    out_dict['plot_layout'] = plot_layout
    return out_dict

####

def extract_above_threshold(df, volcano_results):
    substrate = volcano_results['negfold_pval']['substrate'] + volcano_results['posfold_pval']['substrate']
    print(len(substrate))
    df_substrate = pd.DataFrame.from_dict({'substrate':substrate})
    df2 = df_substrate.join(df.set_index(['substrate']), on=['substrate'])
    return df2

#substrate  control_mean          mean  fold_change   p-value    ctrlCV   treatCV

#####

def KSEA(df, kinase_substrate):
    df.dropna(axis=1, how='all', inplace=True)

    df['subst_position'] = ''
    df['subst_position'] = df.substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    # get name and remove _human
    df['subst_name'] = df.substrate.str.split("(", n=1, expand=True)[0]

    print(df)

    df = df[~(df['fold_change'] == 0)]
    df = df.dropna(axis=0)

    df['Log2substrate_fold_change'] = np.log2(df.fold_change)

    df['Log2substrate_fold_change'] = df['Log2substrate_fold_change'].replace([np.inf, -np.inf], np.nan)

    df = df.dropna(axis=0)

    # unique values
    mean_log2_FC = df["Log2substrate_fold_change"].mean()
    standard_deviation = np.std(df.Log2substrate_fold_change)


    kinase_substrate['substrate'] =  kinase_substrate['substrate'].str.split("_HUMAN", n=1, expand=True)[0]
    kinase_substrate['substrate'] +=  "_HUMAN"


    df = df.join(kinase_substrate[['kinase', 'sub_gene', 'sub_mod_rsd']].set_index(['sub_gene', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'])


    df = df.join(kinase_substrate[['kinase', 'substrate', 'sub_mod_rsd']].set_index(['substrate', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'], rsuffix='_sub_acc')

    df['kinase_final'] = df[['kinase', 'kinase_sub_acc']].fillna('').sum(axis=1)
    # kinase counts
    kinase_count = df['kinase_final'].value_counts()
    kinase_count_df = pd.DataFrame(kinase_count)
    kinase_count_df['sqrt_kinase'] = np.sqrt(kinase_count_df)




    #non identif
    non_identified = kinase_count_df.loc['', 'kinase_final']
    kinase_count_df.drop([''], axis=0, inplace=True)
    kinase_count_df['Kinase'] = kinase_count_df.index
    kinase_count_df.reset_index(drop=True, inplace=True)


    # join this
    tmp_df= pd.DataFrame(df.groupby('kinase_final')['Log2substrate_fold_change'].mean())
    tmp_df['kinase_final'] =  tmp_df.index
    tmp_df.reset_index(drop=True, inplace=True)
    kinase_count_df = kinase_count_df.join(tmp_df.set_index(['kinase_final']), on=['Kinase'], rsuffix='_sub_acc')

    kinase_count_df['KSEA'] = ((kinase_count_df['Log2substrate_fold_change'] - mean_log2_FC) *
                               kinase_count_df['sqrt_kinase']) / standard_deviation


    kinase_count_df['P_value'] = norm.sf(abs(kinase_count_df['KSEA']))


    kinase_count_df = kinase_count_df.sort_values(['KSEA'], ascending=True)


    sig_kinase_count = kinase_count_df[kinase_count_df['P_value'] < 0.05]
    colors = ["rgb(255,242,0)"]*sig_kinase_count[sig_kinase_count['KSEA']<0].shape[0] +  ["rgb(0,34,255)"]*sig_kinase_count[sig_kinase_count['KSEA']>=0].shape[0]
    kinase_count_dict={'kinase': list(sig_kinase_count.loc[:,'Kinase']),
                       'ksea' : list(sig_kinase_count.loc[:,'KSEA']),
                       'p_value': list(sig_kinase_count.loc[:,'P_value']),
                       'colors':colors}

    output = {'score': kinase_count_dict ,'non_identified': non_identified}

    return output
