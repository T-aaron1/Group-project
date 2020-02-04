import requests
from bs4 import BeautifulSoup



def get_uniprot_data(accession_id):
    ''' get data from uniprot api. The data is in xml. This uses bs4 BeautifulSoup to parse it
     for eg: https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession=P31749  OR
     https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&taxid=9606&gene=BRCA1  
     taxid 9606 is for human'''
    text = 'https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&taxid=9606&gene={}'.format(accession_id)
    try:
        requested_xml = requests.get(text)
    except:
        print('Problem with the following accession number: ' + accession_id)
    bs_data = BeautifulSoup(requested_xml.content)
    return bs_data



def get_gralinfo(bs_data, accession_id):
    ''' get general info from bs_data from get_uniprot_data
        uniprot_id, full_prot_name,reverse, chromosome, start_gene_coord, genom_end_coord, sequence,
        added at the beginning the accession_id because the request can be done with the gene or the 
        uniprot id and so on'''
#    sequence =  bs_data.gnentries.gnentry.sequence.text # sequence
    uniprot_id = bs_data.gnentries.gnentry.accession.text # uniprot accession 
#    full_prot_name = bs_data.gnentries.gnentry.protein.recommendedname.fullname.text # full prot name
    reverse = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['reverse_strand'] # rev strand
    chromosome = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['chromosome'] #chromosome
#    start_gene_coord = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['start'] #start
#    genom_end_coord = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['end'] #ends
#    text_list = [accession_id, uniprot_id, full_prot_name,reverse, chromosome, start_gene_coord, genom_end_coord, sequence]
    text_list = [accession_id, uniprot_id, reverse, chromosome]
    text_output = '|'.join(text_list)
    return text_output


# need to extract short name, some do not have short name     
def get_alternative_prot_names(bs_data):
    ''' get alternative prot name from bs_data from get_uniprot_data,
        returns a list'''
    uniprot_id = bs_data.gnentries.gnentry.accession.text
    # get list of all alternative names
    all_alternatives = bs_data.gnentries.gnentry.protein.findAll('alternativename')
    tmp_alternatives = [] # null list
    # append to list
    for alternative in all_alternatives:
        text_list = [uniprot_id, alternative.text]
        text_output = '|'.join(text_list)
        tmp_alternatives.append(text_output)
    return tmp_alternatives
    
#<alternativeName><fullName>Protein kinase B alpha</fullName><shortName>PKB alpha</shortName></alternativeName><alternativeName><fullName>Protein kinase B</fullName><shortName>PKB</shortName>





### bellow are lists
def get_ensembl_geneid(bs_data):
    ''' get ensembl gene id name from bs_data from get_uniprot_data,
        returns a list'''
    uniprot_id = bs_data.gnentries.gnentry.accession.text # uniprot accession 
    genes_ids = []
    genes = bs_data.gnentries.gnentry.findAll('gncoordinate')
    for gene in genes:
        ensembl_gene_id = gene['ensembl_gene_id']
        ensembl_transcript_id = gene['ensembl_transcript_id']
        ensembl_translation_id = gene['ensembl_translation_id']
        print(ensembl_gene_id)
        text_list = [uniprot_id, ensembl_gene_id, ensembl_transcript_id, ensembl_translation_id]
        text_output = '|'.join(text_list)
        genes_ids.append(text_output)

    return genes_ids





## get modified residues

def data_modif_res(bs_data):
    ''' returns a xml list of modif residues, with all information'''
    reference = bs_data.gnentries.gnentry.findAll('gncoordinate')[0]
    features = reference.findAll('feature')
    phospho_list = []
    for entry in features:
        if entry['type'] == 'modified residue': # extract modified residues
             if 'hospho' in entry.findAll('ns2:description')[0].text: # is phospho...something ?
                 phospho_list.append(entry)
    return phospho_list



def modres_gral_info(bs_data):
    ''' General information of Modified residues: position, type, genom_begin, genom_end.
       Columns separated by "|"
       Returns a list. '''
    mod_res_list = data_modif_res(bs_data)
    uniprot_id = bs_data.gnentries.gnentry.accession.text
    tmp_list = []
    
    for mod_res_item in mod_res_list:
        
        residue_position = str(mod_res_item.findAll('ns2:position')[0]['position']) # position
        type_modif = mod_res_item.findAll('ns2:description')[0].text # type of modif
        
        genom_begin = str(mod_res_item.findAll('ns2:begin')[0]['position']) #genom_beg_pos
        genom_end = str(mod_res_item.findAll('ns2:end')[0]['position']) #genom_end_pos
        
        text_list = [uniprot_id, residue_position, type_modif, genom_begin, genom_end]
        text_row = '|'.join(text_list)
        
        tmp_list.append(text_row)
        
    return tmp_list



def modres_references(bs_data):
    ''' References of Modified residues: residue_position, id,type 
       Columns separated by "|"
       Returns a list. '''
    uniprot_id = bs_data.gnentries.gnentry.accession.text
    mod_res_list = data_modif_res(bs_data)
    tmp_main_list = []
    
    for mod_res_item in mod_res_list:
        residue_position = mod_res_item.findAll('ns2:position')[0]['position'] # position
        
        references = mod_res_item.findAll('ns2:dbreference') # position of modif residue
        tmp_ref_list = []
        for reference in references:
            tmp_tup = (reference['id'],reference['type']) # (ref_id,ref_type)
            tmp_ref_list.append(tmp_tup)
        
        for item in tmp_ref_list:
            text_list = [uniprot_id, residue_position, item[0], item[1]]
            text_output = '|'.join(text_list)
            tmp_main_list.append(text_output)
    
    return tmp_main_list

#data_modif_res = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].findAll('feature')
#for entry in data_modif_res:
#    if entry['type'] == 'modified residue': #only mod res
#        if 'hospho' in entry.findAll('ns2:description')[0].text: #description is phosphosomething
#           print(entry.findAll('ns2:description')[0].text) # what is mod. needs to order this
#            print(entry.findAll('ns2:position')[0]['position']) # position of modification
#            print(entry.findAll('ns2:evidence')[0].findAll('ns2:dbreference')[0]['id']) #needs a loop #pubmedid
#            print(entry.findAll('ns2:evidence')[0].findAll('ns2:dbreference')[0]['type'])#needs a loop #pubmed
