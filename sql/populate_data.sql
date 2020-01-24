-- sqlite3 commands

--x families.csv
-- family_abbreviation,residue,family_name
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/families.csv families

SELECT COUNT(*) FROM kinase_info;
SELECT * FROM families;

--x kinase_general_information.csv
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_general_information.csv kinase_info


SELECT COUNT(*) FROM kinase_info;



--x kinase_isoforms_uniprotxml.csv
-- uniprot|tmp_iso_id
.separator "|"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_isoforms_uniprotxml.csv isoforms

SELECT COUNT(*) FROM isoforms;


--x phosphosite_table.csv
-- uniprot_id,residue_position,Modif,type_modif,genom_begin,genom_end
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/phosphosite_table.csv phosphosites


SELECT COUNT(*) FROM phosphosites;



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
.separator "\t"
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_references.csv references_full

SELECT COUNT(*) FROM references_full;


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


--x Kinase_Substrate_Dataset.csv
-- EIF2AK1,HRI,Q9BQI3,eIF2-alpha,1965,P05198,EIF2S1,S52,447635,MILLsELsRRRIRsI,S1, ,X
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/Kinase_Substrate_Dataset.csv  kinase_substrate

SELECT COUNT(*) FROM kinase_substrate;



-- kinase_modified_residues_references_ref_numbers.csv
-- P31749,450,73
.separator ","
.import /homes/dtg30/Desktop/group_proj/venv/src/Group\-project/csv_tables/kinases/kinase_modified_residues_references_ref_numbers.csv  phosphosites_references

SELECT COUNT(*) FROM phosphosites_references;

