import pandas as pd


#substrates_file = "/homes/dtg30/Desktop/Kinase_Substrate_Dataset.csv"
substrates_file = "/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/Kinase_Substrate_Dataset.csv"

data = pd.read_csv(substrates_file)
data.dropna(axis= 1, how='all', inplace = True)


#upload_file = "/homes/dtg30/Desktop/az20.tsv"
upload_file = "/home/daniel/Escritorio/az20.tsv"
upload = pd.read_csv(upload_file, sep = '\t')


upload.dropna(axis= 1, how='all', inplace = True)


upload['subst_position'] = ''
upload['subst_position'] = upload.Substrate.str.split("(", n=1, expand = True)[1].str.replace(")","")
upload['subst_name'] = upload.Substrate.str.split("(", n=1, expand = True)[0]


#get kinase name for substrate / position

# possibly solve this with a join, instead of a for loop, mainly beacause of the time it takes

upload[[]] 
data.columns
upload['kinase'] = ''


upload2 = upload
#upload2 = upload2.join(data[['KINASE','SUBSTRATE','SUB_MOD_RSD']].set_index(['SUBSTRATE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )

upload2 = upload2.join(data[['KINASE','SUB_GENE','SUB_MOD_RSD']].set_index(['SUB_GENE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )

upload2 = upload2.join(data[['KINASE','SUB_GENE','SUB_ACC_ID']].set_index(['SUB_GENE','SUB_ACC_ID']), on= ['subst_name','subst_position'] , rsuffix='_sub_acc')


upload2.columns
upload2['kinase_final'] = upload2[['KINASE','KINASE_sub_acc']].fillna('').sum(axis=1)





            




#add kinase column
upload['kinase'] = upload




#
df = pd.DataFrame({'key1': ['K0', 'K1', 'K2', 'K2', 'K4', 'K5'],
                   'key12':[1,2,3,4,5,6],
                   'A': ['A0', 'A1', 'A2', 'A3', 'A4', 'A5']})
other = pd.DataFrame({'key2': ['K0', 'K2', 'K2'],
                      'key22':[1,9,3],
                      'B': ['B0', 'B1', 'B2']})

df.join(other[['key2','B','key22']].set_index(['key2','key22']), on=['key1','key12'])



dat1 = pd.DataFrame({'dat1': [[9,5,1],[9,5,1],[9,5,1]], 'bla':[3,2,3]})
dat1

# split list, generate rows
dat1.dat1.apply(pd.Series).merge(dat1, left_index=True, right_index = True)\
    .drop(['dat1'], axis=1)\
    .melt(id_vars = ['bla'], value_name='ingg')\
    .drop('variable', axis=1)

df1 = pd.DataFrame({'A': ['a', 'c'], 'B': [np.nan, 'b']})
df2 = pd.DataFrame({'A': [1, 1], 'B': [3, 3]})
take_smaller = lambda s1, s2: s1 if s1.sum() < s2.sum() else s2
df1.combine(df2, take_smaller)
take_smaller = lambda s1, s2: s2 if bool(s1) else s1
df1[['A','B']].combine(df2, take_smaller)

df1[['A','B']].fillna('').sum(axis=1)
