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



SELECT kin.*,basic.chromosome, basic.gene,basic.reverse  FROM kinase_info kin LEFT JOIN basic_info basic ON basic.uniprot_id = kin.uniprot_id WHERE kin.uniprot_id LIKE "P31749";

SELECT subs.sub_acc_id, subs.site_7_aa,subs.sub_mod_rsd,basic.gene AS sub_gene, basic.chromosome,
unip.genom_begin ||':'||unip.genom_end
FROM kinase_substrate subs
LEFT JOIN basic_info basic ON basic.uniprot_id = subs.sub_acc_id
LEFT JOIN uniprot_phosphosites unip ON unip.uniprot_id = subs.sub_acc_id AND unip.residue_position = SUBSTR(subs.sub_mod_rsd,2)
WHERE kin_acc_id LIKE 'P31749'; 

kinase_substrate subs LEFT JOIN basic_info basic ON basic.uniprot_id = subs.sub_acc_id LEFT JOIN uniprot_phosphosites unip ON unip.uniprot_id = subs.sub_acc_id AND unip.residue_position = SUBSTR(subs.sub_mod_rsd,2)

SELECT kin.prot_sequence,basic.chromosome, basic.reverse
FROM kinase_info kin LEFT JOIN basic_info basic ON basic.uniprot_id = kin.uniprot_id

SELECT basick.prot_name AS kinase, basicsub.gene AS sub_gene, subs.sub_mod_rsd, basicsub.prot_name AS substrate
FROM kinase_substrate subs
LEFT JOIN basic_info basick ON subs.kin_acc_id = basick.uniprot_id
LEFT JOIN  basic_info basicsub ON subs.sub_acc_id = basicsub.uniprot_id;

SELECT ncbi.ncbi_id, kin.reverse, subst.genom_begin, subst.genom_end, kin.uniprot_id|| "("||subst.residue_position||")"
  FROM kinase_info kin
  LEFT JOIN ncbi_chrom_id ncbi ON kin.chromosome = ncbi.chr
  LEFT JOIN phosphosites subst ON kin.uniprot_id = subst.uniprot_id
  WHERE kin.chromosome LIKE "x" AND subst.genom_begin NOT NULL;


SELECT basic.reverse, upho.* FROM
uniprot_phosphosites upho
LEFT JOIN basic_info basic ON basic.uniprot_id = upho.uniprot_id
WHERE basic.chromosome LIKE 1 AND upho.genom_begin > 100 AND upho.genom_end <2000000;

