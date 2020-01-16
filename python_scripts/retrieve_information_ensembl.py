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
