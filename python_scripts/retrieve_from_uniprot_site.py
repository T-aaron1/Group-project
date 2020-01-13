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
