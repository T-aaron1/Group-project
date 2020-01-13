#for retrieve1
import retrieve_from_xml_uniprotapi as xmluniprot

uniprot = 'P31749'

data = xmluniprot.get_uniprot_data(uniprot)


grainfo = xmluniprot.get_gralinfo(data)

alt_prot = xmluniprot.get_alternative_prot_names(data)

ensembl_id = xmluniprot.get_ensembl_geneid(data)

modifres_gralinfo = xmluniprot.modres_gral_info(data)

modres_refs = xmluniprot.modres_references(data)


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
