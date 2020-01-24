import sqlite3
import pandas as pd

db = sqlite3.connect('/homes/dtg30/Desktop/group_proj_2/kinase_project.db')
cursor = db.cursor()

cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(cursor.fetchall())

#[('families',), ('sqlite_sequence',), ('kinase_info',), ('isoforms',), ('reactions',), ('references_full',), ('subcell_location_text',),
#('reactions_refences',), ('function_references',), ('kin_function',), ('phosphosites',), ('diseases',), ('subcell_location',),
# ('kinase_substrate',), ('phosphosites_references',)]

cursor.execute("SELECT * from families")
table = pd.read_sql_query("SELECT * from families", db)

table = pd.read_sql_query("SELECT *  FROM families WHERE family_abbreviation LIKE 'ADCK'", db)

#this gives the number of results
pd.read_sql_query("SELECT *  FROM families WHERE family_abbreviation LIKE 'A'", db).shape[0]


# table name
cursor.execute("SELECT * FROM kinase_info")
[item[0] for item in cursor.description]
#['uniprot_id', 'full_prot_name', 'reverse', 'chromosome', 'prot_sequence', 'fasd_name', 'prot_name', 'name_human', 'mass', 'ensembl_gene_id', 'genome_starts', 'genome_ends', 'genome_sequence']

query = "SELECT uniprot_id FROM kinase_info WHERE uniprot_id LIKE '{}'".format('P31749')
pd.read_sql_query(query, db).shape[0]

q_prot_seq = "SELECT prot_sequence FROM kinase_info WHERE uniprot_id LIKE '{}'".format('P31749')
text = pd.read_sql_query(q_prot_seq, db)
text.loc[0,'prot_sequence']
