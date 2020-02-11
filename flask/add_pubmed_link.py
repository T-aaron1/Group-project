import re

def pubmed_link(text):
    '''This function is used to modify the text of some sections in the kinase 
       data (like function and cellular location), which sometimes contain 
       data from Pubmed publication. It adds an <a></a> (link) element so that
       the html text contains hyperlinks to the actual Pubmed reference. The
       function is used inside flask scripts.
    '''
    pmid = re.findall('(?<=PubMed:)[0-9]+', text)
    list_without_pmid = re.split('(?<=PubMed:)[0-9]+', text)
    text_out = list_without_pmid[0]
    for i in range(len(pmid)):
        tmp_link = "<a href='https://www.ncbi.nlm.nih.gov/pubmed/{0}'>{0}</a> ".format(pmid[i])
        text_out += tmp_link
        text_out += list_without_pmid[i+1]
    return text_out
