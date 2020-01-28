# Importing the csv's where data is stored in

import pandas as pd
from pandas import DataFrame

conn = sqlite3.connect('DraftKinaseDB.db')  
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

# Importing disease csv file 
df5 = pd.read_csv('diseases.csv')

# Importing kinase_function csv file 
df6 = pd.read_csv('kinase_function.csv')

# Importing kinase_substrate csv file 
df7 = pd.read_csv('kinase_substrate.csv')

# Importing reactions csv file 
df8 = pd.read_csv('reactions.csv')

# Importing subcellular_location csv file 
df9 = pd.read_csv('subcellular_location.csv')

# Importing subcellular_additional_text csv file 
df10 = pd.read_csv('subcellular_additional_text.csv')
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

# Replace the values from the csv file into the table 'DISEASE'
df5.to_sql('DISEASE', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'KINASE_FUNCTION'
df6.to_sql('KINASE_FUNCTION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'KINASE_SUBSTRATE'
df7.to_sql('KINASE_SUBSTRATE', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'REACTIONS'
df8.to_sql('REACTIONS', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'SUBCELLULAR_LOCATION'
df9.to_sql('SUBCELLULAR_LOCATION', conn, if_exists='replace', index = False)

# Replace the values from the csv file into the table 'SUBCELLULAR_ADDITIONAL_TEXT'
df10.to_sql('SUBCELLULAR_ADDITIONAL_TEXT', conn, if_exists='replace', index = False)
