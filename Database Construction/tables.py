import sqlite3
sqlite3.connect('DraftKinaseDB.db') # Define connection to sqlite3 server

conn = sqlite3.connect('DraftKinaseDB.db')  # You can create a new database by changing the name within the quotes
c = conn.cursor() # The database will be saved in the location where your 'py' (python script) file is saved

# To make database, make the tables first

# Create table - KINASE
c.execute('''CREATE TABLE KINASE
            ([uniprot_id] VARCHAR PRIMARY KEY, [full_prot_name] VARCHAR, [reverse] TEXT, [chromosome] VARCHAR, [sequence] VARCHAR, [Family] TEXT, [name] VARCHAR, [name_human] VARCHAR, [mass] FLOAT, [ensembl_gene_id] VARCHAR, [genome_starts] INTEGER, [genome_ends] INTEGER, [genome_sequence] TEXT)''')

# Create table - PHOSPHOSITE
c.execute('''CREATE TABLE PHOSPHOSITE
             ([phosposite_uniprot_id] VARCHAR PRIMARY KEY, [residue_position] INTEGER, [Modification] TEXT, [Type_Modification] VARCHAR, [Genome_begin] INTEGER, [Genome_end] INTEGER)''')

# create FAMILY table
c.execute('''CREATE TABLE FAMILY
            ([family_abbreviation] VARCHAR PRIMARY KEY, [residue] VARCHAR, [family_name] VARCHAR)''')

# create ISOFORM table
c.execute('''CREATE TABLE ISOFORM
            ([isoform_uniprot] VARCHAR PRIMARY KEY, [tmp_iso_id] VARCHAR)''')

# create DISEASE table
c.execute('''CREATE TABLE DISEASE
            ([uniprot] VARCHAR PRIMARY KEY, [disease_name] VARCHAR, [effect_text] TEXT, [disease_description] TEXT)''')

# create KINASE_FUNCTION table
c.execute('''CREATE TABLE KINASE_FUNCTION
            ([function_uniprot_id] VARCHAR PRIMARY KEY, [protein_fucntion] VARCHAR)''')

# Create table - KINASE_SUBSTRATE
c.execute('''CREATE TABLE KINASE_SUBSTRATE
            ([sub_acc_id] VARCHAR PRIMARY KEY, [gene] VARCHAR, [kinase] VARCHAR, [kin_acc_id] VARCHAR, [substrate] VARCHAR, [sub_gene_id] VARCHAR, [sub_gene] VARCHAR, [sub_mod_rsd] VARCHAR, [site_grp_id] INTEGER, [site_7_aa] VARCHAR, [prot_domain] VARCHAR, [in_vivo_rxn] TEXT, [in_vitro_rxn] TEXT)''')

# create REACTIONS table
c.execute('''CREATE TABLE REACTIONS
            ([reaction_uniprot_id] VARCHAR PRIMARY KEY, [react_id] VARCHAR, [reaction_text] VARCHAR)''')

# create SUBCELLULAR_LOCATION table
c.execute('''CREATE TABLE SUBCELLULAR_LOCATION 
            ([subcellular_uniprot_id] VARCHAR PRIMARY KEY, [subcellular_location] TEXT)''')

# create SUBCELLULAR_ADDITIONAL_TEXT table
c.execute('''CREATE TABLE SUBCELLULAR_ADDITIONAL_TEXT 
            ([subcellular_additional_text_uniprot_id] VARCHAR PRIMARY KEY, [subcellular_aditional_text] TEXT)''')

conn.commit()
