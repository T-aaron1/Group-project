import re

def pubmed_link(text):
    pmid = re.findall('(?<=PubMed:)[0-9]+', text)
    list_without_pmid = re.split('(?<=PubMed:)[0-9]+', text)
    text_out = list_without_pmid[0]
    for i in range(len(pmid)):
        tmp_link = "<a href='https://www.ncbi.nlm.nih.gov/pubmed/{0}'>{0}</a> ".format(pmid[i])
        text_out += tmp_link
        text_out += list_without_pmid[i+1]
    return text_out
