import retrieve_from_xml_uniprotapi as xmluniprot

csvFile =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_list.csv'
csvGralInfo =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_gral_info.csv'
csvAltProt =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_alternative_names.csv'
csvEnsembl =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_ensembl.csv'
csvModResGral =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_modified_residues_gral_info.csv'
csvModResRefs =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_modified_residues_references.csv'

INFILE = open(csvFile, 'r')

OUT_GRALINFO = open(csvGralInfo,'w')
OUT_ALTPROT = open(csvAltProt,'w')
OUT_ENSEMBL = open(csvEnsembl,'w')
OUT_MODRESGRAL = open(csvModResGral,'w')
OUT_MODRESREFS = open(csvModResRefs,'w')

OUT_GRALINFO.write('uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence')
OUT_ALTPROT.write('uniprot_id|name|short')
OUT_ENSEMBL.write('uniprot_id|ensembl_gene_id|ensembl_transcript_id|ensembl_translation_id')
OUT_MODRESGRAL.write('uniprot_id|residue_position|type_modif|genom_begin|genom_end')
OUT_MODRESREFS.write('uniprot_id|residue_position|ref_id|ref_type')




# OUT_GRALINFO = open(,'w')


from

uniprot = 'P31749'

data = xmluniprot.get_uniprot_data(uniprot)

gralinfo = xmluniprot.get_gralinfo(data)
for entry in gralinfo:
    OUT_GRALINFO.

alt_prot = xmluniprot.get_alternative_prot_names(data)
alt_prot

ensembl_id = xmluniprot.get_ensembl_geneid(data)
ensembl_id

modifres_gralinfo = xmluniprot.modres_gral_info(data)
modifres_gralinfo

modres_refs = xmluniprot.modres_references(data)
modres_refs
