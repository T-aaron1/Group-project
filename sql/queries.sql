-- name of all the tables
.tables

-- names of the columns of the tables

.schema diseases
select disease_description from diseases;

SELECT tmp_refs FROM diseases WHERE uniprot LIKE 'P31749';

SELECT uniprot_id, prot_name, name_human FROM kinase_info WHERE uniprot_id LIKE 'P31749' OR prot_name LIKE 'P31749'  OR name_human LIKE 'AKT1_HUMAN';


uniprot_id prot_name name_human
--
