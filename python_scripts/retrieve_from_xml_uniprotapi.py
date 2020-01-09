def get_uniprot_data(accession):
    ''' get data from uniprot api. The data is in xml. This uses bs4 BeautifulSoup to parse it
     for eg: https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession=P31749'''
    text = 'https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession={}'.format(accession)
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


def get_accession(bs_data):
    ''' get accession number from bs_data from get_uniprot_data'''
    return bs_data.gnentries.gnentry.accession.text


def get_sequence(bs_data):
    ''' get sequence from bs_data from get_uniprot_data'''
    return bs_data.gnentries.gnentry.sequence.text

def get_full_prot_name(bs_data):
    ''' get full protein name from bs_data from get_uniprot_data'''
    return bs_data.gnentries.gnentry.protein.recommendedname.fullname.text


def get_alternative_prot_names(bs_data):
    ''' get alternative prot name from bs_data from get_uniprot_data,
        returns a list'''
    alternative_names = [] # null list

    # get list of all alternative names
    all_alternatives = bs_data.gnentries.gnentry.protein.findAll('alternativename')

    # append to list
    for alternative in all_alternatives:
        alternative_names.append(alternative.text)
    return alternative_names
    # should return a list

def get_rev(bs_data):
    ''' get if it is on the reverse strand from bs_data from get_uniprot_data,
        takes the first sequence as reference and extracts info from there
        returns a True/False'''
    reverse = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['reverse_strand']) # rev strand
    return reverse


def get_chrom(bs_data):
    ''' get chromosome from bs_data from get_uniprot_data,
        takes the first sequence as reference and extracts info from there
        returns a number'''
    chromosome = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['chromosome']) #chromosome
    return chromosome


def get_genm_coord_start(bs_data):
    ''' get gene genom coordinate from bs_data from get_uniprot_data,
        takes the first sequence as reference and extracts info from there
        returns a number'''
    start_gene_coord = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['start']) #start
    return start_gene_coord


def get_genm_coord_end(bs_data):
    ''' get gene genom coordinate from bs_data from get_uniprot_data,
        takes the first sequence as reference and extracts info from there
        returns a number'''
    genom_end_coord = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].genomiclocation['end']) #ends
    return genom_end_coord



### bellow are lists

def get_ensembl_geneid(bs_data):
    ''' get ensembl gene id name from bs_data from get_uniprot_data,
        returns a list'''
    genes_ids = []
    genes = bs_data.gnentries.gnentry.findAll('gncoordinate')
    for gene in genes:
        gene_ids.append(gene['ensembl_gene_id'] )
    return genes_ids
    # should return a list


def get_ensembl_transcid(bs_data):
    ''' get ensembl transcript id name from bs_data from get_uniprot_data,
        returns a list'''
    genes = bs_data.gnentries.gnentry.findAll('gncoordinate')
    gene_transcids = []
    for gene in genes:
        gene_transcids.append(gene['ensembl_transcript_id'])
    return gene_transcids
    # should return a list


def get_ensembl_tranlid(bs_data):
    ''' get ensembl translation id from bs_data from get_uniprot_data,
        returns a list'''
    genes = bs_data.gnentries.gnentry.findAll('gncoordinate')
    gene_translids = []
    for gene in genes:
        gene_translids.append(gene['ensembl_translation_id'])
    return gene_translids
    # should return a list


## get modified residues

def data_modif_res(bs_data):
    ''' returns a xml list of modif residues, with all information'''
    reference = bs_data.gnentries.gnentry.findAll('gncoordinate')[0]
    features = reference.findAll('feature')
    phospho_list = []
    for entry in features:
        if entry['type'] == 'modified residue': # extract modified residues
             if 'hospho' in entry.findAll('ns2:description')[0].text: # is phosphosomething ?
                 phospho_list.append(entry)
    return phospho_list



#data_modif_res = bs_data.gnentries.gnentry.findAll('gncoordinate')[0].findAll('feature')
#for entry in data_modif_res:
#    if entry['type'] == 'modified residue': #only mod res
#        if 'hospho' in entry.findAll('ns2:description')[0].text: #description is phosphosomething
#           print(entry.findAll('ns2:description')[0].text) # what is mod. needs to order this
#            print(entry.findAll('ns2:position')[0]['position']) # position of modification
#            print(entry.findAll('ns2:evidence')[0].findAll('ns2:dbreference')[0]['id']) #needs a loop #pubmedid
#            print(entry.findAll('ns2:evidence')[0].findAll('ns2:dbreference')[0]['type'])#needs a loop #pubmed
