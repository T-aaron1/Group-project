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


############################################################################

data = get_uniprot_data('P31749')

uniprot = data.uniprot.entry.findAll('accession')[0].text
data.uniprot.entry.findAll('name').findAll('')
bla = data.findAll('comment')[0]
bla
data.findAll('reaction')[0]['type']
data.findAll('comment')[1].findAll('text')[0].text
[0]['evidence']


###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
### KINASE_UNIPROT_FULL_REF, KINASE_FUNCTION_REFS, KINASE_FUNCTION, REACTIONS_ID NOT CHANGING

file_kinase_function = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_uniprotxml.csv'
file_kinase_function_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_function_refs_uniprotxml.csv'
file_kinase_reactions = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_uniprotxml.csv'
file_kinase_reactions_refs = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_reactions_refs_uniprotxml.csv'
file_kinase_references = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_uniprot_full_references_uniprotxml.csv'
file_kinase_isoforms = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_isoforms_uniprotxml.csv'
file_kinase_diseases = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_diseases_uniprotxml.csv'
file_kinase_subcell_loc = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_uniprotxml.csv'
file_kinase_subcell_loc_text = '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_subcell_loc_text_uniprotxml.csv'



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
csvFile =  '/homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/kinase_list_uniprotxml.csv'
INFILE = open(csvFile, 'r')
counter = 1

print(INFILE.readline())

error_get_prot_function = []
error_get_diseases= []
error_get_kinase_isoforms= []
error_get_reactions= []
error_get_full_references= []
error_get_subcell_loc= []

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





OUTFILE_FUNCTION_REFS.close()
OUTFILE_FUNCTION.close()
OUTFILE_DISEASES.close()
OUTFILE_kinase_isoforms.close()
OUTFILE_kinase_reactions.close()
OUTFILE_kinase_reactions_refs.close()
OUTFILE_kinase_references.close()
OUTFILE_kinase_subcell_loc.close()
OUTFILE_kinase_subcell_loc_text.close()



### muchos errores, especialmente en get_full_references,
###   algunas veces en: get_subcell_loc, get_prot_function, get_reactions
### Q8TD19:
# 218:  Q8TD19
# !! issue: Q8TD19: get_reactions
# !! issue: Q8TD19: get_full_references
# !! issue: Q8TD19: get_subcell_loc
# 204:  Q13523
# !! issue: Q13523: get_prot_function
# !! issue: Q13523: get_reactions
# !! issue: Q13523: get_full_references
# !! issue: Q13523: get_subcell_loc
# 202:  Q9UQ07
# !! issue: Q9UQ07: get_full_references
# !! issue: Q9UQ07: get_subcell_loc
# 292:  P15056
# !! issue: P15056: get_full_references


###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# functions



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





# @@@@@@
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
            print(function_text_output)
            OUTFILE_FUNCTION.write(function_text_output)
            prot_function_refs = comment.findAll('text')[0]['evidence']
            prot_function_refs_list += prot_function_refs.split()
    for item in prot_function_refs_list:
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        text_out_list = [uniprot, item]
        text_output = '|'.join(text_out_list) + '\n'
        OUTFILE_FUNCTION_REFS.write(text_output)

#get_prot_function(data, OUTFILE_FUNCTION, OUTFILE_FUNCTION_REFS)

#OUTFILE_FUNCTION_REFS.close()
#OUTFILE_FUNCTION.close()

# @@@@@@@@@@@@
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
                effect_text = comment.findAll('text')[0].text
            except:
                effect_text = ''
            try:
                disease_description = disease_info.findAll('description')[0].text
            except:
                disease_description = ''
            text_out_list = [uniprot, tmp_refs, disease_id, disease_name, effect_text, disease_description]
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
    reaction_text_refs_list = []
    reaction_list = data.findAll('reaction')
    reaction_refs_ids = []
    for item in reaction_list :
        uniprot = data.uniprot.entry.findAll('accession')[0].text
        id_react_refs = uniprot + '_' +str(counter)
        reaction_refs_ids.append(id_react_refs)
        reaction_text = item.findAll('text')[0].text
        tmp_refs_list = item['evidence']
        reaction_text_refs_list.append(tmp_refs_list.split())
        text_out_list = [id_react_refs, reaction_text]
        text_output = '|'.join(text_out_list) + '\n'
        OUTFILE.write(text_output)
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
    for reference in daaata: #data.findAll('reference'):
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
            for item in subcell_list:
                info = item.findAll('location')[0]
                try:
                    subcell_refs = info['evidence']
                except:
                    subcell_refs = ''
                subcell_location = info.text
                text_out_list = [uniprot, subcell_location, subcell_refs]
                text_output = '|'.join(text_out_list) +'\n'
                OUTFILE.write(text_output)
            subcell_aditional_text = comment.findAll('text')[0].text
            subcell_aditional_text_refs = comment.findAll('text')[0]['evidence']
            text_out_list = [uniprot, subcell_aditional_text, subcell_aditional_text_refs]
            text_output = '|'.join(text_out_list) + '\n'
            print(text_output)
            OUTFILE_TEXT.write(text_output)

# get_subcell_loc(data, OUTFILE_kinase_subcell_loc, OUTFILE_kinase_subcell_loc_text)
#
# OUTFILE_kinase_subcell_loc.close()
# OUTFILE_kinase_subcell_loc_text.close()
