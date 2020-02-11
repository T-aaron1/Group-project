-- name of all the tables
.tables

-- names of the columns of the tables

.schema diseases
select disease_description from diseases;

SELECT tmp_refs FROM diseases WHERE uniprot LIKE 'P31749';

SELECT uniprot_id, prot_name, name_human FROM kinase_info WHERE uniprot_id LIKE 'P31749' OR prot_name LIKE 'P31749'  OR name_human LIKE 'AKT1_HUMAN';


uniprot_id prot_name name_human
--


-- inhibitors
inhibitors_gral_info --x
inhibitors_kin_family --x
inhibitors_pdbid   -- x
inhibitors_synonims -- x
inhibitors_targets -- x

SELECT targets FROM inhibitors_targets WHERE inn_name LIKE "Afuresertib"


SELECT targets FROM inhibitors_targets WHERE inn_name LIKE "Afuresertib" JOIN kinase_info.uniprot_id;


SELECT inh.targets, kin.uniprot_id FROM inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets WHERE inn_name LIKE "Afuresertib" ;

.schema

ncbi_chrom_id ncbi_id

chromosome  kinase_info


SELECT ncbi.ncbi_id, kin.reverse, subst.genom_begin, subst.genom_end, kin.uniprot_id|| "("||subst.residue_position||")"
  FROM kinase_info kin
  LEFT JOIN ncbi_chrom_id ncbi ON kin.chromosome = ncbi.chr
  LEFT JOIN phosphosites subst ON kin.uniprot_id = subst.uniprot_id
  WHERE kin.chromosome LIKE "x" AND subst.genom_begin NOT NULL;
