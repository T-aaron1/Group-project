### TO JOIN THE TABLES USING A FOREIGN KEY

# Joining KINASE table to PHOSPHOSITE TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN PHOSPHOSITE ON PHOSPHOSITE.phosposite_uniprot_id = KINASE.uniprot_id ;''') 
#print(cur.fetchall())

# Joining KINASE table to FAMILY TABLE
cur.execute('''
SELECT Family FROM KINASE INNER JOIN FAMILY ON FAMILY.family_abbreviation = KINASE.Family ;''') 

# Joining KINASE table to ISOFORM TABLE
cur.execute('''
SELECT uniprot_id FROM KINASE INNER JOIN ISOFORM ON ISOFORM.isoform_uniprot = KINASE.uniprot_id ;''') 

# Joining KINASE table to DISEASE TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN DISEASE ON DISEASE.uniprot = KINASE.uniprot_id ;''') 

# Joining KINASE table to KINASE_FUNCTION TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN KINASE_FUNCTION ON KINASE_FUNCTION.function_uniprot_id = KINASE.uniprot_id ;''') 

# Joining KINASE table to KINASE_SUBSTRATE TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN KINASE_SUBSTRATE ON KINASE_SUBSTRATE.sub_acc_id = KINASE.uniprot_id ;''')

# Joining KINASE table to REACTIONS TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN REACTIONS ON REACTIONS.reaction_uniprot_id = KINASE.uniprot_id ;''') 

# Joining KINASE table to SUBCELLULAR_LOCATION TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN SUBCELLULAR_LOCATION ON SUBCELLULAR_LOCATION.subcellular_uniprot_id = KINASE.uniprot_id ;''') 

# Joining KINASE table to SUBCELLULAR_ADDITIONAL_TEXT TABLE
cur.execute('''SELECT uniprot_id FROM KINASE INNER JOIN SUBCELLULAR_ADDITIONAL_TEXT ON SUBCELLULAR_ADDITIONAL_TEXT.subcellular_additional_text_uniprot_id = KINASE.uniprot_id ;''') 

