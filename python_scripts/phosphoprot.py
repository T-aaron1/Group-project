import pandas as pd

substrates_file = "/homes/dtg30/Desktop/Kinase_Substrate_Dataset.csv"

data = pd.read_csv(substrates_file)
data.dropna(axis= 1, how='all', inplace = True)
data.columns


upload_file = "/homes/dtg30/Desktop/az20.tsv"
upload = pd.read_csv(upload_file, sep = '\t')

upload.dropna(axis= 1, how='all', inplace = True)


upload['subst_position'] = upload.Substrate.str.split("(", n=1, expand = True)[1].str.replace(")","")
upload['subst_name'] = upload.Substrate.str.split("(", n=1, expand = True)[0]


SUB_ORGANISM
SUB_GENE_ID
SUB_MOD_RSD


#get kinase name for substrate / position

for subs in zip(upload['subst_name'], upload['subst_position']):
    if (subs[0] in list(data['SUBSTRATE'])):
       print(list (data[(data['SUBSTRATE'] == subs[0]) & (data['SUB_MOD_RSD'] == subs[1])]['KINASE']))
    elif (subs[0] in list(data['SUBSTRATE'])):
       print(list(data[(data['SUBSTRATE'] == subs[0]) & (data['SUB_MOD_RSD'] == subs[1])]['KINASE']))


upload['kinase'] = ''

for i in range(len(upload)):
    if (upload['subst_name'][i] in list(data['SUBSTRATE'])):
        tmp_df = data[data['SUBSTRATE'] == upload['subst_name'][i]]
        kinase = list(tmp_df[tmp_df['SUB_MOD_RSD'] == upload['subst_position'][i]]['KINASE'])
        if kinase:
            print('substrate')
            print(kinase)
            upload['kinase'][i] = kinase
            print(upload['kinase'][i])
    elif (upload['subst_name'][i] in list(data['SUB_GENE'])):
        tmp_df = data[data['SUB_GENE'] == upload['subst_name'][i]]
        kinase = list(tmp_df[tmp_df['SUB_MOD_RSD'] == upload['subst_position'][i]]['KINASE'])
        if kinase:
            print('subgene')
            print(kinase)
            upload['kinase'][i] = kinase
            print(upload['kinase'][i])






#add kinase column
upload['kinase'] = upload
