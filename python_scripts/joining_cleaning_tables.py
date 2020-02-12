'''
scripts used to join and clean some csv files
short scripts to handle particular things
'''

import pandas as pd


#### to remove columns and split files to create general_basic_info table
import pandas as pd
# uniprot_id,full_prot_name,reverse,chromosome,sequence,Family,name,name_human,mass,
#ensembl_gene_id,genome_starts,genome_ends,genome_sequence
INFILE_KINASES = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinase_general_information.csv'
INFILE_TARGET_gral = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets/kinase_gral_info.csv'
PHOSPLUS = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/Kinase_Substrate_Dataset.csv'
OUTPUT_FILE_BASICINFO = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/general_information_basic.csv'
OUTPUT_FILE_KININFO = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinase_general_kin_information.csv' 

phosplus = pd.read_csv(PHOSPLUS)
infilekin = pd.read_csv(INFILE_KINASES)
infilekin.columns
phosplus.columns

target_chrrev = pd.read_csv(INFILE_TARGET_gral, sep='|')

basic_info = infilekin[['uniprot_id','reverse','chromosome', 'name']]
gral_kininfo = infilekin.drop(columns = ['reverse','chromosome', 'name'])

phosplus.drop(columns = ['sub_gene_id', 'site_grp_id'], inplace=True) 
phosplus.drop(columns = ['prot_domain', 'in_vivo_rxn','in_vitro_rxn'], inplace=True)
# whole columns
#['gene', 'kinase', 'kin_acc_id',
#'substrate', 'sub_acc_id', 'sub_gene',
#       'sub_mod_rsd', 'site_7_aa'],

#wanted 
#['gene', 'kinase', 'kin_acc_id',
#'substrate', 'sub_acc_id', 'sub_gene',
#: uniprot_id, gene, phosphositeplusname  ... need to add reverse, chromosome from target_chrrev
phosplus.columns

df1 = phosplus[['gene', 'kinase', 'kin_acc_id']]
df2 = phosplus[['substrate', 'sub_acc_id', 'sub_gene']]
df1.rename(columns = {'kin_acc_id':'uniprot_id', 'kinase':'phosphplus_name'}, inplace=True)
df2.rename(columns = {'substrate':'phosphplus_name', 'sub_acc_id':'uniprot_id', 'sub_gene':'gene'}, inplace=True)

df_combined = df1.append(df2)
df_combined.drop_duplicates(inplace=True)
df_combined.columns
target_chrrev.columns
df_combined_chrev = df_combined.join(target_chrrev.set_index('uniprot_id'), on=['uniprot_id'])
df_combined_chrev.columns
# combine with basic_info
basic_info.columns
basic_info.rename(columns={'name':'gene'}, inplace=True)

lista_ya = list(basic_info['uniprot_id'])

basic_info.columns
basic_info_phplus = basic_info.join(df_combined_chrev[['phosphplus_name', 'uniprot_id']].set_index('uniprot_id'), on =['uniprot_id'])

df_comb_chrev_notused = df_combined_chrev[~df_combined_chrev.uniprot_id.isin(lista_ya)]

basic_info_combined = basic_info_phplus.append(df_comb_chrev_notused)
basic_info_combined.columns
basic_info_combined.shape[0]
basic_info_combined.rename(columns={ 'phosphplus_name':'prot_name'}, inplace=True)
#think this is what I wanted.. prot_name not used for the queries, but needed for the analysis as substrate
basic_info_combined.columns


#wanted2:
#'kin_acc_id',
#'sub_acc_id',
#       'sub_mod_rsd', 'site_7_aa'],
# need to add genomic location... but with a join in sql, not in the table
phosplus.columns
phosphosites = phosplus[['kin_acc_id','sub_acc_id','sub_mod_rsd', 'site_7_aa']]
phosphosites.columns
'kin_acc_id', 'sub_acc_id', 'sub_mod_rsd', 'site_7_aa'

gral_kininfo.columns
phosphosites.columns
basic_info_combined.columns

path = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/'
gral_kininfo.to_csv(path+'kin_gral_info_updated.csv', index=False)
phosphosites.to_csv(path+'kin_phosphosites_updated.csv', index=False)
basic_info_combined.to_csv(path+'kin_basic_info_updated.csv', index=False)


### to include phosphosites in phosphosites tables
# remove the ones that are already there
import pandas as pd
import sqlite3

INFILE_CSV = "/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets/kinase_modified_residues_gral_info_2.csv"

OUTPUT_FILE = "/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets/kinase_modified_residues_gral_info_filtered2.csv"

data = pd.read_csv(INFILE_CSV, sep='|')

DATABASE = "/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinase_project.db"
db = sqlite3.connect(DATABASE)
db_results = pd.read_sql_query("SELECT * FROM uniprot_phosphosites", db)
db.close()

uniprot_ids = db_results.uniprot_id.unique()

data.columns
data_filtered = data[~data.uniprot_id.isin(uniprot_ids)]
data_filtered.to_csv(OUTPUT_FILE, index=False)





#####################
#
import pandas as pd
INFILE_CSV = "~/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc.csv"
data = pd.read_csv(INFILE_CSV)
data.shape[0]
data.drop_duplicates().shape[0]
data2 = data.drop_duplicates()
data2.to_csv("~/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc_nodupl.csv", index=False)

#####################
# inhibitors

OUTPUT_FILE = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_table_unique_information.csv'

inhib_gralinfo_csv = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Cleaned_inhibitors_table_information.csv'
inhib_chemform_csv = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitors_inchikey_smiles_unique.csv'


inhib_gralinfo = pd.read_csv(inhib_gralinfo_csv)
inhib_chemform = pd.read_csv(inhib_chemform_csv)

inhib_gralinfo.columns
inhib_chemform.columns

# INN_Name
inhib_gralinfo = inhib_gralinfo.join(inhib_chemform.set_index(['INN_Name']), on = ['INN_Name'])
inhib_gralinfo.to_csv(OUTPUT_FILE, index = False)

########## cleaning the kin_substrate_dataset
data_csv =  '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv'

OUTPUT_FILE = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinases/gene_names.csv'

data = pd.read_csv(data_csv)

# 'KINASE', 'KIN_ACC_ID', 'GENE'
genes_networkin = data.GENE[data.Source == 'NetworKIN'].drop_duplicates()
genes_networkin.shape

genes_networkin.to_csv(OUTPUT_FILE)












OUTPUT_FILE = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_general_information.csv'

######################
## ensembl id cleaning

ensembl_id_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_ensembl.csv'
#uniprot_id|ensembl_gene_id|ensembl_transcript_id|ensembl_translation_id

ensembl_id = pd.read_csv(ensembl_id_csv, sep= '|')

ensembl_id.columns


ee2 = ensembl_id[['uniprot_id', 'ensembl_gene_id']].drop_duplicates()

dupl_list = list(ee2.uniprot_id[ee2['uniprot_id'].duplicated()].drop_duplicates())


for i in range(len(dupl_list)):
    ee2[(ee2['uniprot_id'] == dupl_list[i])]

dupl_list_ensembl_main = ['ENSG00000117020','ENSG00000185974','ENSG00000072062','ENSG00000117676','ENSG00000139908','ENSG00000213923','ENSG00000134058','ENSG00000176444','ENSG00000105204','ENSG00000181085','ENSG00000104814','ENSG00000129465','ENSG00000204580','ENSG00000146904','ENSG00000174292','ENSG00000198870','ENSG00000214102','ENSG00000173137']

for i in range(len(dupl_list)):
    index_drop = ee2[(ee2['uniprot_id'] == dupl_list[i]) &  (ee2['ensembl_gene_id'] != dupl_list_ensembl_main[i])].index
    ee2.drop(index_drop,inplace = True)

ee2.shape
ee2.columns
#['uniprot_id', 'ensembl_gene_id']
### ee2 contains all the main sequences for the genes.

gene_info_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_gene_information.csv'
#uniprot|ensembl_id|genome_starts|genome_ends|genome_sequence

#read and drop duplicates
gene_info = pd.read_csv(gene_info_csv, sep = '|')
gene_info = gene_info.drop_duplicates()
gene_info.shape
gene_info.columns
# ['uniprot', 'ensembl_id', 'genome_starts', 'genome_ends', 'genome_sequence']

gene_information = ee2.join(gene_info.set_index(['uniprot', 'ensembl_id']), on = ['uniprot_id', 'ensembl_gene_id'])

# gene information of all 504 genes is in gene_information Dataframe
######################

gral_info_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_gral_info.csv'
#uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence


family_mass_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_list.csv'
#Family,name,name_human,uniprot,mass

function_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_uniprotxml.csv'
#uniprot|prot_function

disease_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_diseases_uniprotxml.csv'
#uniprot|tmp_refs|disease_id|disease_name|effect_text|disease_description



gral_info = pd.read_csv(gral_info_csv, sep='|')
family_mass = pd.read_csv(family_mass_csv, sep = ',')
function = pd.read_csv(function_csv, sep = '|')

final = gral_info.drop(['start_gene_coord', 'genom_end_coord'], axis =1)
final.shape
final = final.join(family_mass.set_index(['uniprot']), on = ['uniprot_id'] )
final.shape

final.columns


#check the thing of the genes, it would be better to get the list and retrieve the info again..... join with gene info
final = final.join(gene_information.set_index(['uniprot_id']), on = ['uniprot_id'] , how='left')
final.shape



final.to_csv(OUTPUT_FILE, index= False)

# different table
final  = final.join(function.set_index(['uniprot']), on = ['uniprot_id'] , how='')
final.shape


#########################
# clean kin_subcell_loc.csv

OUTPUT_FILE = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc_not_dupl.csv'
kin_subcell_loc_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc.csv'

kin_subcell_loc = pd.read_csv(kin_subcell_loc_csv)

kin_subcell_loc.drop_duplicates(inplace=True)
kin_subcell_loc.to_csv(OUTPUT_FILE)


#########################
# clean kin_subcell_loc_text.csv

OUTPUT_FILE = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc_text_non_duplicates.csv'
kin_subcell_loc_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kin_subcell_loc_text.csv'

kin_subcell_loc = pd.read_csv(kin_subcell_loc_csv)
kin_subcell_loc.shape
kin_subcell_loc.drop_duplicates(inplace=True)
kin_subcell_loc.shape
kin_subcell_loc.to_csv(OUTPUT_FILE)



#############################
# add ref numbers kinase_modified_residues_references.csv
mod_res_refs_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_modified_residues_references.csv'
OUTPUT_FILE = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_modified_residues_references_ref_numbers.csv'
fullrefs_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_references.csv'

mod_res_refs = pd.read_csv(mod_res_refs_csv)
fullrefs = pd.read_csv(fullrefs_csv, sep='\t')

tmp_df = fullrefs[['uniprot', 'reference_id', 'pubmedid']]
tmp_df = tmp_df[~np.isnan(tmp_df.pubmedid)  ]

tmp_df = tmp_df.astype({'pubmedid':int})
tmp_df = tmp_df.astype({'pubmedid':str})

mod_res_refs = mod_res_refs.astype({'ref_id': str})


mod_res_refs = mod_res_refs.join(tmp_df.set_index('pubmedid'), on = ['ref_id'])

mod_res_refs = tmp_df[~np.isnan(mod_res_refs.reference_id)]
mod_res_refs = mod_res_refs.astype({'reference_id':int})
mod_res_refs.drop(['uniprot','ref_type','ref_id'],axis=1, inplace=True)
mod_res_refs.columns
mod_res_refs.to_csv(OUTPUT_FILE, index = False)
