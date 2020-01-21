import pandas as pd

gral_info_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_gral_info.csv'
#uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence

gene_info_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_gene_information.csv'
#uniprot|ensembl_id|genome_starts|genome_ends|genome_sequence

family_mass_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_list.csv'
#Family,name,name_human,uniprot,mass

function_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_uniprotxml.csv'
#uniprot|prot_function

disease_csv = '~/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_diseases_uniprotxml.csv'
#uniprot|tmp_refs|disease_id|disease_name|effect_text|disease_description



gral_info = pd.read_csv(gral_info_csv, sep='|')
gene_info = pd.read_csv(gene_info_csv, sep = '|')
gene_info = gene_info.drop_duplicates()
family_mass = pd.read_csv(family_mass_csv, sep = ',')
function = pd.read_csv(function_csv, sep = '|')


gene_info = gene_info.drop_duplicates()

final = gral_info.drop(['start_gene_coord', 'genom_end_coord'], axis =1)
final.shape
final = final.join(family_mass.set_index(['uniprot']), on = ['uniprot_id'] )
final.shape



# different table
final  = final.join(function.set_index(['uniprot']), on = ['uniprot_id'] , how='')
final.shape

#check the thing of the genes, it would be better to get the list and retrieve the info again
final = final.join(gene_info.set_index(['uniprot']), on = ['uniprot_id'] , how='left')
final.shape
