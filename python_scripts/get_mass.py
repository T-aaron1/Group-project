# python3
'''short script to get kinase mass from Uniprot API'''
import requests

INFILE= open('/homes/dtg30/Desktop/group_proj/venv/src/Group-project/kinase_list.csv','r')
OUTFILE = open('/homes/dtg30/Desktop/group_proj/venv/src/Group-project/python_scripts/kinase_list_mass.csv','w')

header = INFILE.readline()

header = header + ',mass' + '\n'
OUTFILE.write(header)

for line in INFILE:
    line = line.rstrip()
    list= line.split(',')
    name = list[2]
    text = 'https://www.uniprot.org/uniprot/?query={}&columns=mass&format=tab'.format(name)
    data = requests.get(text).text
    mass = data.split('\n')[1]
    mass = mass.replace(',','')
    new_line = line + ',' + mass +'\n'
    OUTFILE.write(new_line)

OUTFILE.close()
INFILE.close()
    
    









