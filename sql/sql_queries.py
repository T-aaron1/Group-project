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
