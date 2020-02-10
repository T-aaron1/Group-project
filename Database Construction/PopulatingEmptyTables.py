# Importing the csv's where data is stored in

import pandas as pd
from pandas import DataFrame

conn = sqlite3.connect('final_version_kinase.db')  
c = conn.cursor()

# Importing kinase csv file 
df1 = pd.read_csv('kinase.csv')
#print(read_kinase)

# Importing phosphosite csv file 
df2 = pd.read_csv('phosphosite.csv')

# Importing family csv file 
df3 = pd.read_csv('families.csv')

# Importing isoform csv file 
df4 = pd.read_csv('kinase_isoforms.csv')

# Importing isoform additional information csv file 
df5 = pd.read_csv('isoform_additional_info.csv')

# Importing disease csv file 
df6 = pd.read_csv('diseases.csv')

# Importing kinase_function csv file 
df7 = pd.read_csv('kinase_function.csv')

# Importing kinase_substrate csv file 
df8 = pd.read_csv('kinase_substrate.csv')

# Importing reactions csv file 
df9 = pd.read_csv('reactions.csv')

# Importing subcellular_location csv file 
df10 = pd.read_csv('subcellular_location.csv')

# Importing subcellular_additional_text csv file 
df11 = pd.read_csv('subcellular_additional_text.csv')

# Importing inhibitor_general_information  csv file 
df12 = pd.read_csv('inhibitors_gral_info.csv')

# Importing inhibitor_kinase_family csv file 
df13 = pd.read_csv('inhibitors_kin_family.csv')

# Importing inhibitor_pdbid csv file 
df14 = pd.read_csv('inhibitors_pdbid.csv')

# Importing inhibitor_synonyms csv file 
df15 = pd.read_csv('inhibitors_synonyms.csv')

# Importing inhibitor_targets csv file 
df16 = pd.read_csv('inhibitors_targets.csv')
_______________________________________________________________________________________________________________________________


# Populating the empty tables 

# Insert the values from the csv file into the table 'KINASE'
df1.to_sql('KINASE', conn, if_exists='append', index = False)

# Replace the values from the csv file into the table 'PHOSPHOSITE'
df2.to_sql('PHOSPHOSITE', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'FAMILY'
df3.to_sql('FAMILY', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'ISOFORM'
df4.to_sql('ISOFORM', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'ISOFORM_ADDITIONAL_INFORMATION'
df5.to_sql('ISOFORM_ADDITIONAL_INFORMATION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'DISEASE'
df6.to_sql('DISEASE', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'KINASE_FUNCTION'
df7.to_sql('KINASE_FUNCTION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'KINASE_SUBSTRATE'
df8.to_sql('KINASE_SUBSTRATE', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'REACTIONS'
df9.to_sql('REACTIONS', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'SUBCELLULAR_LOCATION'
df10.to_sql('SUBCELLULAR_LOCATION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'SUBCELLULAR_ADDITIONAL_TEXT'
df11.to_sql('SUBCELLULAR_ADDITIONAL_TEXT', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'INHIBITOR_GENERAL_INFORMATION'
df12.to_sql('INHIBITOR_GENERAL_INFORMATION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'INHIBITOR_KINASE_FAMILY'
df13.to_sql('INHIBITOR_KINASE_FAMILY', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'INHIBITOR_PDBID'
df14.to_sql('INHIBITOR_PDBID', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'INHIBITOR_SYNONYMS'
df15.to_sql('INHIBITOR_SYNONYMS', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'INHIBITOR_TARGETS'
df16.to_sql('INHIBITOR_TARGETS', conn, if_exists='replace', index = False)
