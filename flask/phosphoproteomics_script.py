import os
import pandas as pd
import numpy as np
from scipy.stats import norm



def file_handle(file_path):
    ''' reads phosphoproteomics data, returns a dictionary with two entries: substrate-position and
        another entry that contains the foldchange information for the inhibitors. Each inhibitor dataframe is
        inside a dictionary with the name of the inhibitor. This output is to be used as input of other functions
        for phosphoproteomics analysis.
    '''
    df = pd.DataFrame(pd.read_csv(file_path, sep='\t'))
    os.remove(file_path)
    df.columns = map(str.lower, df.columns)
    #drop rows and columns with no value
    df.dropna(axis=1, how='all', inplace=True)
    df.dropna(axis=0, how='all', inplace=True)
    #split substrate column into name and position 
    df['subst_position'] = df.substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    df['subst_name'] = df.substrate.str.split("(", n=1, expand=True)[0]
    substrate = df[['substrate','subst_name','subst_position']] # only those columns
    #extracting unique inhibitor name 
    inhibitors_list = []
    inhibitors_dict={}
    for name in df.columns:
        if ('fold_change' in name):
                inhibitor = name.rsplit('_')[0]
                inhibitor = inhibitor.lower().rstrip()
                inhibitors_list.append(inhibitor)
    for i in range(len(inhibitors_list)):
        tmp_inhib_df=df.filter(regex= inhibitors_list[i])
        inhibitor = inhibitors_list[i] + '_'
        tmp_inhib_df.columns = [inh.replace(inhibitor,'') for inh in tmp_inhib_df.columns]
        if 'treatcv' in tmp_inhib_df.columns :
            inhibitors_dict[inhibitors_list[i]]= tmp_inhib_df[['fold_change', 'p-value','treatcv']] # only those columns
        else:
            inhibitors_dict[inhibitors_list[i]]= tmp_inhib_df[['fold_change', 'p-value']] # only those columns

    output ={'substrate': substrate, 'inhibitors_dict' :inhibitors_dict }
    return output


#df instead of file_path...
def volcano(df, p_val_threshold, fold_threshold, cv_threshold): #!! change this, put in dif file
    pval_threshold = float(p_val_threshold)
    fold_threshold = float(fold_threshold)
    cv_threshold = float(cv_threshold)
    
    if 'treatcv' in df.columns :
        df = df[df['treatcv'] < cv_threshold ]

        
    list_nulls = list(df[df.fold_change.isnull()].substrate)
    df = df[(df.fold_change.notnull()) &  (df.fold_change != 0)]
    df['log_foldchange'] = np.log2(df['fold_change'])
    df['log_pval'] = -np.log10(df['p-value'])

    df =df[~(np.isinf(np.abs(df.log_foldchange)))& ~(np.isinf(np.abs(df.log_pval))) & ~np.isnan(df.log_foldchange) & ~np.isnan(df.log_pval)]


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

    df1 = df.join(kinase_substrate[['kinase', 'sub_gene', 'sub_mod_rsd']].set_index(['sub_gene', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'])

    df2 = df.join(kinase_substrate[['kinase', 'substrate', 'sub_mod_rsd']].set_index(['substrate', 'sub_mod_rsd']),
                 on=['subst_name', 'subst_position'])

    df = df1.append(df2)
    df.drop_duplicates(inplace=True)

    df['kinase'] = df[['kinase']].fillna('')
    df.drop_duplicates(inplace=True)
    
    kinase_count = df['kinase'].value_counts()
    kinase_count_df = pd.DataFrame(kinase_count)
    kinase_count_df['sqrt_kinase'] = np.sqrt(kinase_count_df)

    print("df3=============")
    print(kinase_count_df)

    
    non_identified = kinase_count_df.loc['', 'kinase']
    kinase_count_df.drop([''], axis=0, inplace=True)
    kinase_count_df['Kinase'] = kinase_count_df.index
    kinase_count_df.reset_index(drop=True, inplace=True)

    # join this
    tmp_df= pd.DataFrame(df.groupby('kinase')['Log2substrate_fold_change'].mean())
    tmp_df['kinase'] =  tmp_df.index
    tmp_df.reset_index(drop=True, inplace=True)
    kinase_count_df = kinase_count_df.join(tmp_df.set_index(['kinase']), on=['Kinase'], rsuffix='_sub_acc')

    kinase_count_df['KSEA'] = ((kinase_count_df['Log2substrate_fold_change'] - mean_log2_FC) *
                               kinase_count_df['sqrt_kinase']) / standard_deviation

    kinase_count_df['P_value'] = norm.sf(abs(kinase_count_df['KSEA']))
    kinase_count_df = kinase_count_df.sort_values(['KSEA'], ascending=True)
    kinase_count_df2 =  kinase_count_df.round(2)
    sig_kinase_count = kinase_count_df2[kinase_count_df2['P_value'] < 0.05]

    colors = ["rgb(255,242,0)"]*sig_kinase_count[sig_kinase_count['KSEA']<0].shape[0] + \
        ["rgb(0,34,255)"]*sig_kinase_count[sig_kinase_count['KSEA']>=0].shape[0]

    kinase_count_dict={'kinase': list(sig_kinase_count.loc[:,'Kinase']),
                       'ksea' : list(sig_kinase_count.loc[:,'KSEA']),
                       'p_value': list(sig_kinase_count.loc[:,'P_value']),
                       'colors':colors,
                       }

    output = {'score': kinase_count_dict ,'non_identified': non_identified}

    return output
