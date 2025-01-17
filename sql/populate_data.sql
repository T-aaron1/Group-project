-- sqlite3 commands

--x families.csv
-- family_abbreviation,residue,family_name
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/families.csv families

SELECT COUNT(*) FROM kinase_info;
SELECT * FROM families;




--x kinase_isoforms_uniprotxml.csv
-- uniprot|tmp_iso_id
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_isoforms_uniprotxml.csv isoforms

SELECT COUNT(*) FROM isoforms;



--x kinase_general_information.csv
-- new
-- kin_gral_info_updated.csv
-- uniprot_id,full_prot_name,sequence,family,name_human,mass,ensembl_gene_id,genome_starts,genome_ends,genome_sequence
.separator ","
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/kinases/kin_gral_info_updated.csv kinase_info


SELECT COUNT(*) FROM kinase_info;
-- new
-- kin_phosphosites_updated.csv 
-- prev:: phosphosite_table.csv
-- uniprot_id,residue_position,Modif,type_modif,genom_begin,genom_end
.separator ","
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/kinases/phosphosite_table.csv  uniprot_phosphosites

-- /home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/targets/kinase_modified_residues_gral_info_filtered.csv
-- uniprot_id,residue_position,modif,type,genom_begin,genom_end
.separator "|"
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/kinases/targets/kinase_modified_residues_gral_info_filtered2.csv  uniprot_phosphosites

SELECT COUNT(*) FROM uniprot_phosphosites;



-- diseases.csv
-- uniprot_id,residue_position,Modif,type_modif,genom_begin,genom_end
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/diseases.csv diseases

SELECT COUNT(*) FROM diseases;



--x kinase_function_uniprotxml.csv
-- uniprot|prot_function
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_function_uniprotxml.csv kin_function

SELECT COUNT(*) FROM kin_function;


--x kinase_function_refs_uniprotxml.csv
-- P31749|11
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_function_refs_uniprotxml.csv function_references

SELECT COUNT(*) FROM function_references;


--x kinase_reactions_uniprotxml.csv
-- P31749|1|ATP + L-seryl-[protein] = ADP + H(+) + O-phospho-L-seryl-[protein]
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_reactions_uniprotxml.csv reactions

SELECT COUNT(*) FROM reactions;


--x kinase_reactions_refs_uniprotxml.csv
-- 37|P31749|1
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_reactions_refs_uniprotxml.csv reactions_refences

SELECT COUNT(*) FROM reactions_refences;


--x kinase_references.csv
-- P31749	1	journal article	1991	Proc. Natl
--  .separator "\t"
-- .import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_references.csv references_full

-- SELECT COUNT(*) FROM references_full;


--x kin_subcell_loc_text.csv
-- P31749,41.0,
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kin_subcell_loc_text.csv subcell_location_text

SELECT COUNT(*) FROM subcell_location_text;




--x kin_subcell_loc.csv
-- P31749,Cytoplasm,11
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kin_subcell_loc.csv  subcell_location

SELECT COUNT(*) FROM subcell_location;


-- new
-- kin_phosphosites_updated.csv
--previous: Kinase_Substrate_Dataset.csv
--previous: EIF2AK1,HRI,Q9BQI3,eIF2-alpha,1965,P05198,EIF2S1,S52,447635,MILLsELsRRRIRsI,S1, ,X
.separator ","
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/kinases/kin_phosphosites_updated.csv  kinase_substrate 

SELECT COUNT(*) FROM kinase_substrate;

-- new
-- kin_basic_info_updated.csv
--chromosome,gene,prot_name,reverse,uniprot_id
.separator ","
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/kinases/kin_basic_info_updated.csv basic_info

SELECT COUNT(*) FROM basic_info;

-- kinase_modified_residues_references_ref_numbers.csv
-- P31749,450,73
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_modified_residues_references_ref_numbers.csv  phosphosites_references

SELECT COUNT(*) FROM phosphosites_references;







-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/isoforms/kinase_gral_info.csv
-- uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence
.mode csv
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/isoforms/kinase_gral_info.csv  isoforms_gral_info

SELECT COUNT(*) FROM phosphosites_references;


-- -- csv_tables/ncbi_chromosomes
-- /Escritorio/uk/group_project/venv/src/Group-project/csv_tables/ncbi_chromosomes
-- chr,ncbi_id
.separator ','
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group\-project/csv_tables/ncbi_chromosomes ncbi_chrom_id




--------------------------------------------------
--------------------------------------------------
-- inhibitors

-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_table_unique_information.csv
-- inn_name,phase,mw,image_url,canonical_smiles,inchikey
.separator ','
.import /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_table_unique_information.csv  inhibitors_gral_info


-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_synonims_table.csv
-- inn_name,synonyms
.separator ','
.import /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_synonims_table.csv inhibitors_synonims



-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitor_kinase_families.csv
-- inn_name,kinase_familie
.separator ','
.import /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitor_kinase_families.csv inhibitors_kin_family



-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitors_pdbID.csv
-- inn_name,pdbid
.separator ','
.import /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitors_pdbID.csv inhibitors_pdbid



-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_target.csv
-- inn_name,targets
.separator ','
.import /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_target.csv inhibitors_targets

--/home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinase_alternative_names.csv
-- uniprot_id|name|short
.mode csv
.separator "|"
.import /home/daniel/Escritorio/uk/group_project/venv/src/Group-project/csv_tables/kinases/kinase_alternative_names.csv  kinase_alternative_names
