def KSEA(df, kinase_substrate):

    df.dropna(axis= 1, how='all', inplace = True)


    df['subst_position'] = ''
    df['subst_position'] = df.Substrate.str.split("(", n=1, expand = True)[1].str.replace(")","")
    df['subst_name'] = df.Substrate.str.split("(", n=1, expand = True)[0]

    df = df.dropna(axis=0)


    df = df.join(kinase_substrate[['KINASE','SUB_GENE','SUB_MOD_RSD']].set_index(['SUB_GENE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )

    df= df.join(kinase_substrate[['KINASE','SUB_ACC_ID','SUB_MOD_RSD']].set_index(['SUB_ACC_ID','SUB_MOD_RSD']), on= ['subst_name','subst_position'] , rsuffix='_sub_acc')

    df['kinase_final'] = df[['KINASE','KINASE_sub_acc']].fillna('').sum(axis=1)
    df = df[df.fold_change != 'inf']
    df['Log2Substrate_fold_change'] = np.log2(df.fold_change)
    mean_log2_FC =df["Log2Substrate_fold_change"].mean()
    kinase_count =df['kinase_final'].value_counts()
    kinase_count_df = pd.DataFrame(kinase_count)
    kinase_count_df['sqrt_kinase'] = np.sqrt(kinase_count_df)
    df.groupby('kinase_final')['Log2Substrate_fold_change'].mean()
    kinase_count_df['Standard_deviation']= np.std(df.Log2Substrate_fold_change)
    kinase_count_df['Log2Substrate_fold_change'] =df.groupby('kinase_final')['Log2Substrate_fold_change'].mean()
    kinase_count_df['mean_log2(FC)'] = df['mean_log2(FC)']=df["Log2Substrate_fold_change"].mean()
    kinase_count_df['KSEA'] =  (kinase_count_df['Log2Substrate_fold_change'] - mean_log2_FC *
                            kinase_count_df['sqrt_kinase'])/ kinase_count_df['Standard_deviation']

    kinase_count_df['P_value']= norm.sf(abs(kinase_count_df['KSEA']))*2
    non_identified= kinase_count_df.loc['','kinase_final']
    kinase_count_df.drop([''], axis=0, inplace=True )
    kinase_count_df=kinase_count_df.sort_values(['P_value'], ascending = True)
    output= {'score':kinase_count_df, 'non_identified': non_identified}
    return output
