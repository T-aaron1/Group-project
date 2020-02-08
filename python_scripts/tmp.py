####
# tests join vs append in panda's data frames
import pandas as pd

df1 = df = pd.DataFrame({"a":[1, 2, 3, 4], 
                         "b":[5, 6, 7, 8]}) 
  
# Creating the Second Dataframe using dictionary 
df2 = pd.DataFrame({"a":[1, 2, 3], 
                    "b":[5, 6, 7]}) 

df3 = pd.DataFrame({"a":[3,2,1], 
                    "b":[7,6,5]}) 


# Print df1 
print(df1, "\n") 



# Print df2 
df2

ddf = df1.append(df2)
ddf.drop_duplicates()

ddf2 = df1.append(df3)
ddf2.drop_duplicates()

ddf3 = df1.append(df2, ignore_index = True)
ddf3.drop_duplicates()

ddf.drop_duplicates()
ddf2.drop_duplicates()
ddf3.drop_duplicates()




# Importing pandas as pd 
import pandas as pd 

# Creating the first Dataframe using dictionary 
df1 = pd.DataFrame({"a":[1, 2, 3, 4], 
					"b":[5, 6, 7, 8]}) 

# Creating the Second Dataframe using dictionary 
df2 = pd.DataFrame({"a":[1, 2, 3], 
					"c":[1, 5, 4]}) 

# for appending df2 at the end of df1 
df.append(df2, ignore_index = True) 



#######################################################
# unique targets

INCSVFILE =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv'
import pandas as pd

data = pd.read_csv(INCSVFILE)
data.columns

data_networkin = data[data['source'] == 'NetworKIN']
data.shape
data_networkin.shape
targets = data_networkin.sub_gene.drop_duplicates()
targets.shape
#'kinase', 'kin_acc_id', 'gene', 'kin_organism', 'substrate','sub_gene_id', 'sub_acc_id',
#'sub_gene', 'sub_organism', 'sub_mod_rsd', 'site_grp_id', 'site_7_aa', 'networkin_score', 'source

targets.to_csv("/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/targets2/substrates_targets_networkin.csv", index=False)

#######################################################




TEXTO = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

salida = []
nr = 50
ninr = 10
counter = 10
for i in range(0,len(TEXTO), nr):
    subtexto = TEXTO[i:i+nr]
    tmp_list =[]
    for j in range(0,len(subtexto), ninr):
        tmp_list.append(subtexto[j:j+ninr])
    text_out = ' '.join(tmp_list)
    l_index = []
    for i in range(int(nr/ninr)):
        l_index.append(str(counter + (ninr * i)).rjust(10))
    salida.append(' '.join(l_index))
    salida.append(text_out)
    counter += nr









#    -------------------------------
import pandas as pd

tsv_file = '~/Desktop/conrad.tsv'

data = pd.read_csv(tsv_file, sep = '\t')

data.dropna(axis= 1, how='all', inplace = True)

inhibitor = 'AZ20' + '_'
col_names = data.columns
for name in col_names:
    if inhibitor in name:
       data.rename(columns = {name: name.replace(inhibitor,'')}, inplace = True)

















#for retrieve1
import retrieve_from_xml_uniprotapi as xmluniprot

uniprot = 'P31749'

data =  get_uniprot_data(uniprot)


grainfo =  get_gralinfo(data)

alt_prot =  get_alternative_prot_names(data)

ensembl_id =  get_ensembl_geneid(data)

modifres_gralinfo =  modres_gral_info(data)

modres_refs =  modres_references(data)


# to retrieve from uniprot data site but in xml
import requests
from bs4 import BeautifulSoup
uniprot = 'P31749'
data = get_uniprot_data(uniprot)
text = 'https://www.uniprot.org/uniprot/{}.xml'.format(uniprot)
requested_xml = requests.get(text)
bs_data = BeautifulSoup(requested_xml.content)
if bs_data.uniprot.entry.findAll('comment')[0]['type'] == 'function':

bs_data.uniprot.entry.findAll('comment')[0].findAll('text')[0]['evidence']
bs_data.uniprot.entry.findAll('comment')[0].findAll('text')[0].text




###################################
# sequence divide
sequence = "MEHIQGAWKTISNGFGFKDAVFDGSSCISPTIVQQFGYQRRASDDGKLTDPSKTSNTIRVFLPNKQRTVVNVRNGMSLHDCLMKALKVRGLQPECCAVFRLLHEHKGKKARLDWNTDAASLIGEELQVDFLDHVPLTTHNFARKTFLKLAFCDICQKFLLNGFRCQTCGYKFHEHCSTKVPTMCVDWSNIRQLLLFPNSTIGDSGVPALPSLTMRRMRESVSRMPVSSQHRYSTPHAFTFNTSSPSSEGSLSQRQRSTSTPNVHMVSTTLPVDSRMIEDAIRSHSESASPSALSSSPNNLSPTGWSQPKTPVPAQRERAPVSGTQEKNKIRPRGQRDSSYYWEIEASEVMLSTRIGSGSFGTVYKGKWHGDVAVKILKVVDPTPEQFQAFRNEVAVLRKTRHVNILLFMGYMTKDNLAIVTQWCEGSSLYKHLHVQETKFQMFQLIDIARQTAQGMDYLHAKNIIHRDMKSNNIFLHEGLTVKIGDFGLATVKSRWSGSQQVEQPTGSVLWMAPEVIRMQDNNPFSFQSDVYSYGIVLYELMTGELPYSHINNRDQIIFMVGRGYASPDLSKLYKNCPKAMKRLVADCVKKVKEERPLFPQILSSIELLQHSLPKINRSASEPSLHRAAHTEDINACTLTTSPRLPVF"

seq_out = []
for i in range(0,len(sequence),10):
    seq_out.append(sequence[i:i+10])

cada = 6
seq_out_final = []
for i in range(0,len(seq_out), cada):
    seq_out_final.append(seq_out[i:i+cada])

seq_out_final

'{:9} YES votes  {:2.2%}'.format(yes_votes, percentage)

seq_by_10s = range(0,len(sequence), 10)

'{:-10}{:-10}{:-10}{:-10}{:-10}'.format('1','2','3','4','5')

counter = 10
for entry in seq_out_final:
    str(counter).rjust(10)+str(counter+10).rjust(10)+str(counter+20).rjust(10)+str(counter+30).rjust(10)+str(counter+40).rjust(10)
    print(' '.join(entry))
    counter +=50




##########################
# ensembl

import requests

direction = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/'
file_csv = 'kinase_ensembl.csv'
file_gene_sequence = 'kinase_gene_information.csv'

INFILE = open(direction + file_csv, 'r')
INFILE.readline()
#uniprot_id|ensembl_gene_id|ensembl_transcript_id|ensembl_translation_id\n
OUTFILE_GENE = open(direction + file_gene_sequence, 'w')

OUTFILE_GENE.write('uniprot|ensembl_id|genome_starts|genome_ends|genome_sequence\n')

separator = '|'
ensembl_id_list  = []

for line in INFILE :
    line = line.strip().split('|')
    ensembl_id = line[1]
    if ensembl_id not in ensembl_id_list:
        ddict = {}
        ddict[line[0]] = line[1]
        ensembl_id_list.append(ddict)


for entry in ensembl_id_list:
    uniprot = [*entry][0]
    ensembl_id = entry[uniprot]
    print(ensembl_id)
    server = "https://rest.ensembl.org/sequence/id/{}".format(ensembl_id)
    fasta = requests.get(server, headers={ "Content-Type" : "text/x-fasta"})
    header = fasta.text.split('\n')[0]
    genome_starts = header.split(':')[3]
    genome_ends = header.split(':')[4]
    genome_sequence_request = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    genome_sequence = genome_sequence_request.text.rstrip()
    text_out_list = [uniprot,ensembl_id, genome_starts, genome_ends, genome_sequence]
    text_output = '|'.join(text_out_list) + '\n'
    OUTFILE_GENE.write(text_output)

OUTFILE_GENE.close()
































############################################################################

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###  KINASE_FUNCTION_REFS, KINASE_FUNCTION, REACTIONS_ID NOT CHANGING

# files to create
file_kinase_function = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_uniprotxml.csv'
file_kinase_function_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_refs_uniprotxml.csv'
file_kinase_reactions = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_uniprotxml.csv'
file_kinase_reactions_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_refs_uniprotxml.csv'
file_kinase_references = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_uniprot_full_references_uniprotxml.csv'
file_kinase_isoforms = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_isoforms_uniprotxml.csv'
file_kinase_diseases = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_diseases_uniprotxml.csv'
file_kinase_subcell_loc = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_uniprotxml.csv'
file_kinase_subcell_loc_text = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_text_uniprotxml.csv'

# open all files
OUTFILE_FUNCTION = open(file_kinase_function, 'w')
OUTFILE_FUNCTION_REFS = open(file_kinase_function_refs, 'w')
OUTFILE_DISEASES = open(file_kinase_diseases, 'w')
OUTFILE_kinase_isoforms = open(file_kinase_isoforms, 'w')
OUTFILE_kinase_reactions = open(file_kinase_reactions, 'w')
OUTFILE_kinase_reactions_refs = open(file_kinase_reactions_refs, 'w')
OUTFILE_kinase_references = open(file_kinase_references, 'w')
OUTFILE_kinase_subcell_loc = open(file_kinase_subcell_loc, 'w')
OUTFILE_kinase_subcell_loc_text = open(file_kinase_subcell_loc_text, 'w')


## test
csvFile =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_list.csv'
INFILE = open(csvFile, 'r')
counter = 1

print(INFILE.readline())

# list of errors
error_get_prot_function = []
error_get_diseases= []
error_get_kinase_isoforms= []
error_get_reactions= []
error_get_full_references= []
error_get_subcell_loc= []

# write the headers
OUTFILE_FUNCTION.write('uniprot|prot_function' +'\n')
OUTFILE_FUNCTION_REFS.write('uniprot|item' + '\n')
OUTFILE_DISEASES.write('uniprot|tmp_refs|disease_id|disease_name|eheaderect_text|disease_description' + '\n')
OUTFILE_kinase_isoforms.write('uniprot|tmp_iso_id' + '\n')
OUTFILE_kinase_reactions.write('id_react_refs|reaction_text' + '\n')
OUTFILE_kinase_reactions_refs.write('ref_item|ref_id' + '\n')
OUTFILE_kinase_references.write('uniprot|reference_id|pub_type|pub_date|pub_name|pub_vol|pub_pages_first|pub_pages_last|pub_title|aut_text|pub_dbref_text' + '\n')
OUTFILE_kinase_subcell_loc.write('uniprot|subcell_location|subcell_refs' + '\n')
OUTFILE_kinase_subcell_loc_text.write('uniprot|subcell_aditional_text|subcell_aditional_text_refs' + '\n')

# run all the functions defined bellow to
#  get and extract the data
# try/except to handle non-existent information
for line in INFILE:
    line = line.rstrip()
    tmp_list = line.split(',')
    uniprot = tmp_list[3]
    text = str(counter)+ ':  '+ uniprot
    print(text)
    ## get data
    data = get_uniprot_data(uniprot)
    ## --  gral info
    try:
        get_prot_function(data, OUTFILE_FUNCTION, OUTFILE_FUNCTION_REFS)
    except:
        error_get_prot_function.append(uniprot)
        print("!! issue: " + uniprot + ': ' + 'get_prot_function')
    try:
        get_diseases(data, OUTFILE = OUTFILE_DISEASES)
    except:
        print("!! issue: " + uniprot + ': ' + 'get_diseases')
        error_get_diseases.append(uniprot)
    try:
        get_kinase_isoforms(data, OUTFILE = OUTFILE_kinase_isoforms)
    except:
        print("!! issue: " + uniprot + ': ' + 'get_kinase_isoforms')
        error_get_kinase_isoforms.append(uniprot)
    try:
        get_reactions(data, OUTFILE_kinase_reactions, OUTFILE_kinase_reactions_refs)
    except:
        print("!! issue: " + uniprot + ': ' + 'get_reactions')
        error_get_reactions.append(uniprot)
    try:
        get_full_references(data, OUTFILE_kinase_references)
    except:
        print("!! issue: " + uniprot + ': ' + 'get_full_references')
        error_get_full_references.append(uniprot)
    try:
        get_subcell_loc(data, OUTFILE_kinase_subcell_loc, OUTFILE_kinase_subcell_loc_text)
    except:
        print("!! issue: " + uniprot + ': ' + 'get_subcell_loc')
        error_get_subcell_loc.append(uniprot)
    print("---------------------------")
    counter += 1


# close files
OUTFILE_FUNCTION_REFS.close()
OUTFILE_FUNCTION.close()
OUTFILE_DISEASES.close()
OUTFILE_kinase_isoforms.close()
OUTFILE_kinase_reactions.close()
OUTFILE_kinase_reactions_refs.close()
OUTFILE_kinase_references.close()
OUTFILE_kinase_subcell_loc.close()
OUTFILE_kinase_subcell_loc_text.close()




###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# functions
# comments in all functions to point output_file, run the function, close_output_file

def get_uniprot_data(accession):
    ''' get data from uniprot api. The data is in xml. This uses bs4 BeautifulSoup to parse it
     for eg: https://www.uniprot.org/uniprot/P31749.xml'''
    text = 'https://www.uniprot.org/uniprot/{}.xml'.format(accession)
    try:
        import requests
    except:
        'requests package not installed. Run "pip install -r requirements"'
    try:
        from bs4 import BeautifulSoup
    except:
        print(' bs4 package not installed. Run "pip install -r requirements"')
    try:
        requested_xml = requests.get(text)
    except:
        print('Problem with the following accession number: ' + accession)
    bs_data = BeautifulSoup(requested_xml.content)
    return bs_data

# unique values per prot
###########################################
## function

#file_kinase_function = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_uniprotxml.csv'
#file_kinase_function_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_refs_uniprotxml.csv'

#OUTFILE_FUNCTION = open(file_kinase_function, 'w')
#OUTFILE_FUNCTION_REFS = open(file_kinase_function_refs, 'w')

def get_prot_function(bs_data, OUTFILE_FUNCTION, OUTFILE_FUNCTION_REFS):
    prot_function_refs_list = []
    for comment in data.findAll('comment') :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        if comment['type'] == 'function':
            prot_function = comment.findAll('text')[0].text
            text_out_list = [uniprot, prot_function]
            function_text_output = '|'.join(text_out_list) + '\n'
            OUTFILE_FUNCTION.write(function_text_output)
            try:
                prot_function_refs = comment.findAll('text')[0]['evidence']
            except:
                prot_function_refs = ''
            prot_function_refs_list += prot_function_refs.split()
    for item in prot_function_refs_list:
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        text_out_list = [uniprot, item]
        text_output = '|'.join(text_out_list) + '\n'
        OUTFILE_FUNCTION_REFS.write(text_output)

#get_prot_function(data, OUTFILE_FUNCTION, OUTFILE_FUNCTION_REFS)

#OUTFILE_FUNCTION_REFS.close()
#OUTFILE_FUNCTION.close()

# several
###########################################
## diseases

#file_kinase_diseases = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_diseases_uniprotxml.csv'

#OUTFILE_DISEASES = open(file_kinase_diseases, 'w')

def get_diseases(bs_data, OUTFILE):
    diseases_list = []
    for comment in data.findAll('comment') :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        if comment['type'] == 'disease':
            try:
                tmp_refs = comment['evidence'].replace(' ','/')
            except:
                tmp_refs = ''
            try:
                disease_info = comment.findAll('disease')[0]
            except:
                disease_info = ''
            try:
                disease_id = disease_info['id']
            except:
                disease_id = ''
            try:
                disease_name = disease_info.findAll('name').text
            except:
                disease_name = ''
            try:
                eheaderect_text = comment.findAll('text')[0].text
            except:
                eheaderect_text = ''
            try:
                disease_description = disease_info.findAll('description')[0].text
            except:
                disease_description = ''
            text_out_list = [uniprot, tmp_refs, disease_id, disease_name, eheaderect_text, disease_description]
            text_output = '|'.join(text_out_list) + '\n'
            OUTFILE.write(text_output)

#get_diseases(data, OUTFILE = OUTFILE_DISEASES)

#OUTFILE_DISEASES.close()



###########################################
## isoforms


#file_kinase_isoforms = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_isoforms_uniprotxml.csv'

# OUTFILE_kinase_isoforms = open(file_kinase_isoforms, 'w')

def get_kinase_isoforms(bs_data, OUTFILE):
    for comment in data.findAll('comment') :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        if comment['type'] == 'alternative products':
            tmp_isoforms = comment.findAll('isoform')
            for iso in tmp_isoforms :
                tmp_iso_id = iso.findAll('id')[0].text
                text_out_list = [uniprot, tmp_iso_id]
                text_output = '|'.join(text_out_list)  + '\n'
                OUTFILE.write(text_output)

# get_kinase_isoforms(data, OUTFILE = OUTFILE_kinase_isoforms)
#
# OUTFILE_kinase_isoforms.close()



###########################################
## catalytic activity

# file_kinase_reactions = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_uniprotxml.csv'
# file_kinase_reactions_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_refs_uniprotxml.csv'
#
# OUTFILE_kinase_reactions = open(file_kinase_reactions, 'w')
# OUTFILE_kinase_reactions_refs = open(file_kinase_reactions_refs, 'w')


def get_reactions(bs_data, OUTFILE, OUTFILE_REFS):
    counter = 1
    reaction_text_refs_list = []
    reaction_list = data.findAll('reaction')
    reaction_refs_ids = []
    for item in reaction_list :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        id_react_refs = uniprot + '_' +str(counter)
        reaction_refs_ids.append(id_react_refs)
        reaction_text = item.findAll('text')[0].text
        try:
            tmp_refs_list = item['evidence']
        except:
            tmp_refs_list = ''
        reaction_text_refs_list.append(tmp_refs_list.split())
        text_out_list = [id_react_refs, reaction_text]
        text_output = '|'.join(text_out_list) + '\n'
        OUTFILE.write(text_output)
        counter += 1
    for i in range(len(reaction_refs_ids)):
        ref_id = reaction_refs_ids[i]
        refs_list = reaction_text_refs_list[i]
        for ref_item in refs_list:
            text_out = ref_item + '|' + ref_id + '\n'
            OUTFILE_REFS.write(text_out)


# get_reactions(data, OUTFILE_kinase_reactions, OUTFILE_kinase_reactions_refs)
#

# OUTFILE_kinase_reactions.close()
# OUTFILE_kinase_reactions_refs.close()

###########################################
## references

# file_kinase_references = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_uniprot_full_references_uniprotxml.csv'
#
# OUTFILE_kinase_references = open(file_kinase_references, 'w')

def get_full_references(bs_data, OUTFILE):
    for reference in data.findAll('reference'):
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        reference_id = reference['key']
        ref_data = reference.findAll('citation')[0]
        try:
            pub_type = ref_data['type']
        except:
            pub_type = ''
        try:
            pub_date = ref_data['date']
        except:
            pub_date = ''
        try:
            pub_name = ref_data['name']
        except:
            pub_name = ''
        try:
            pub_vol = ref_data['volume']
        except:
            pub_vol = ''
        try:
            pub_pages_first = ref_data['first']
        except:
            pub_pages_first = ''
        try:
            pub_pages_last = ref_data['last']
        except:
            pub_pages_last = ''
        try:
            pub_title = ref_data.findAll('title')[0].text
        except:
            pub_title = ''
        aut_list = []
        try:
            tmp_author_names = ref_data.findAll('authorlist')[0].findAll('person')
            for author in tmp_author_names:
                aut_list.append(author['name'])
        except:
            pass
        pub_dbref = []
        try:
            tmp_dbrefs = ref_data.findAll('dbreference')
            for dbref in tmp_dbrefs:
                ref_tuple_text = dbref['type']+ ','+ dbref['id']
                pub_dbref.append(ref_tuple_text)
        except:
            pass
        aut_text = '/'.join(aut_list)
        pub_dbref_text = '/'.join(pub_dbref)
        text_out_list = [uniprot,reference_id,pub_type,pub_date,pub_name,pub_vol,pub_pages_first,pub_pages_last,pub_title,aut_text,pub_dbref_text ]
        text_output = '|'.join(text_out_list) + '\n'
        OUTFILE.write(text_output)


# get_full_references(data, OUTFILE_kinase_references)

# OUTFILE_kinase_references.close()


###########################################
## subcellular location

# file_kinase_subcell_loc = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_uniprotxml.csv'
# file_kinase_subcell_loc_text = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_text_uniprotxml.csv'

# OUTFILE_kinase_subcell_loc = open(file_kinase_subcell_loc, 'w')
# OUTFILE_kinase_subcell_loc_text = open(file_kinase_subcell_loc_text, 'w')

def get_subcell_loc(bs_data, OUTFILE, OUTFILE_TEXT):
    for comment in data.findAll('comment') :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        if comment['type'] == 'subcellular location':
            subcell_list = comment.findAll('subcellularlocation')
#            print(subcell_list)
            for item in subcell_list:
                info_list = item.findAll('location')
                for info in info_list:
                    try:
                        subcell_refs = info['evidence']
                    except:
                        subcell_refs = ''
                    subcell_location = info.text
                    text_out_list = [uniprot, subcell_location, subcell_refs]
                    text_output = '|'.join(text_out_list) +'\n'
#                    print(text_output)
                OUTFILE.write(text_output)
            try:
                subcell_aditional_text = comment.findAll('text')[0].text
            except:
                subcell_aditional_text = ''
            try:
                subcell_aditional_text_refs = comment.findAll('text')[0]['evidence']
            except:
                subcell_aditional_text_refs = ''
            text_out_list = [uniprot, subcell_aditional_text, subcell_aditional_text_refs]
            text_output = '|'.join(text_out_list) + '\n'
#            print(text_output)
            OUTFILE_TEXT.write(text_output)

# get_subcell_loc(data, OUTFILE_kinase_subcell_loc, OUTFILE_kinase_subcell_loc_text)
#
# OUTFILE_kinase_subcell_loc.close()
# OUTFILE_kinase_subcell_loc_text.close()
