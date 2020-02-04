import retrieve_from_xml_uniprotapi as xmluniprot

#/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinase_list.csv
#/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets/substrates_list.csv
#/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets2/substrates_targets_networkin.csv
# /home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinases/
###   gene_names.csv


csvFile =  '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinases/gene_names.csv'


path = '/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinases/'
csvGralInfo = path + 'kinase_gral_info.csv'
csvAltProt =  path + 'kinase_alternative_names.csv'
csvEnsembl =  path + 'kinase_ensembl.csv'
csvModResGral =  path + 'kinase_modified_residues_gral_info.csv'
csvModResRefs =  path + 'kinase_modified_residues_references.csv'

INFILE = open(csvFile, 'r')

OUT_GRALINFO = open(csvGralInfo,'w')
OUT_ALTPROT = open(csvAltProt,'w')
OUT_ENSEMBL = open(csvEnsembl,'w')
OUT_MODRESGRAL = open(csvModResGral,'w')
OUT_MODRESREFS = open(csvModResRefs,'w')


# writing headers on output files
OUT_GRALINFO.write('id_in|uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence\n')
OUT_ALTPROT.write('uniprot_id|name|short\n')
OUT_ENSEMBL.write('uniprot_id|ensembl_gene_id|ensembl_transcript_id|ensembl_translation_id\n')
OUT_MODRESGRAL.write('uniprot_id|residue_position|type_modif|genom_begin|genom_end\n')
OUT_MODRESREFS.write('uniprot_id|residue_position|ref_id|ref_type\n')


# read kinase list

# read headers
print(INFILE.readline())
separator = ','
column_wanted = 0

# we want uniprot  [3]
#Family,name,name_human,uniprot,mass
#AGC_Ser/Thr,AKT1,AKT1_HUMAN,P31749,55.686

# errors lists
## --  gral info
error_list_get_gralinfo = []
## -- alternative names
error_get_alternative_prot_names = []
## -- ensembl ids 
error_list_get_ensembl_geneid = []
## -- modif residue information
error_modres_references = []

counter = 1

for line in INFILE:
    line = line.rstrip()
    tmp_list = line.split(separator)
    acc_id_in = tmp_list[column_wanted]
    print('------------------------------------------------------')
    print('----------------')
    print(counter)
    print(acc_id_in)
    print('----------------')
    ## get data
    data = xmluniprot.get_uniprot_data(acc_id_in)
    ## --  gral info
    try:
        gralinfo = xmluniprot.get_gralinfo(data, acc_id_in)
        out =  gralinfo + '\n'
        print(out)
        OUT_GRALINFO.write(out) # write 
    except:
        error_list_get_gralinfo.append(acc_id_in)
        print('!! some info not found: get_gralinfo: '+ acc_id_in)
    ## -- alternative names
    try:
        alt_prot = xmluniprot.get_alternative_prot_names(data)
        if (len(alt_prot) > 1):
            for entry in alt_prot:
                out = entry + '\n'
                OUT_ALTPROT.write(out) #write 
        else:
            out = ''.join(alt_prot) + '\n'
            OUT_ALTPROT.write(out) #write
    except:
        error_get_alternative_prot_names.append(acc_id_in)
        print('!! some info not found: get_alternative_prot_names: '+ acc_id_in)
    ## -- ensembl ids 
    try: 
        ensembl_id = xmluniprot.get_ensembl_geneid(data)
        if(len(ensembl_id) > 1):
            for entry in ensembl_id:
                out = entry + '\n'
                OUT_ENSEMBL.write(out) #write
        else:
            out = ''.join(ensembl_id) + '\n'
            OUT_ENSEMBL.write(out) # write
    except:
        error_list_get_ensembl_geneid.append(acc_id_in)
        print('!! some info not found: get_ensembl_geneid:'+ acc_id_in)
    ## -- modif residue information
    try:
        modifres_gralinfo = xmluniprot.modres_gral_info(data)
        for entry in modifres_gralinfo:
            out = entry + '\n'    
            OUT_MODRESGRAL.write(out) #write
        ## modified residues references
        modres_refs = xmluniprot.modres_references(data)
        if(len(modres_refs) > 1 ): 
            for entry in modres_refs:
                out = entry + '\n'
                OUT_MODRESREFS.write(out) #write
        else: 
            out = ''.join(modres_refs) + '\n'
            OUT_MODRESREFS.write(out) #write
    except:
        error_modres_references.append(acc_id_in)
        print('!! some info not found: modres_references: '+ acc_id_in)
    print('------------------------------------------------------')
    counter += 1


OUT_GRALINFO.close()
OUT_ALTPROT.close()
OUT_ENSEMBL.close()
OUT_MODRESGRAL.close()
OUT_MODRESREFS.close()

INFILE.close()



# errors lists
## --  gral info
error_list_get_gralinfo 
## -- alternative names
error_get_alternative_prot_names 
## -- ensembl ids 
error_list_get_ensembl_geneid 
## -- modif residue information
error_modres_references
