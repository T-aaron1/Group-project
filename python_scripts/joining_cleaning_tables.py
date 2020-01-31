import pandas as pd



########## cleaning the kin_substrate_dataset
















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
