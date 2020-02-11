-- get schema: https://dbdiagram.io/d
-- possibly not all tables were populated
-- run again the comands to populate the tables


-- sqlite3 commands

-- families.csv
CREATE TABLE families (
	   family_abbreviation TEXT PRIMARY KEY,
residue TEXT,
family_name TEXT
);



-- kinase_general_information.csv
CREATE TABLE kinase_info (
    uniprot_id TEXT PRIMARY KEY,
    full_prot_name TEXT,
    reverse TEXT,
    chromosome TEXT,
    prot_sequence TEXT,
    family TEXT, -- segundo cambio ñ
    gene TEXT,  -- primer cambio ñ
    name_human TEXT,
    mass INTEGER,
    ensembl_gene_id TEXT,
    genome_starts INTEGER,
    genome_ends INTEGER,
    genome_sequence TEXT,
    FOREIGN KEY(family) REFERENCES families(family_abbreviation)
);




-- kinase_isoforms_uniprotxml.csv
CREATE TABLE isoforms (
uniprot TEXT,
isoform TEXT PRIMARY KEY,
FOREIGN KEY(uniprot)  REFERENCES kinase_info(uniprot_id)
);


-- phosphosite_table.csv
CREATE TABLE phosphosites (
uniprot_id TEXT,
residue_position INTEGER,
modif TEXT,
type_modif TEXT,
genom_begin INTEGER,
genom_end  INTEGER,
FOREIGN KEY(uniprot_id) REFERENCES kinase_info(uniprot_id)
);


-- uniprot|tmp_refs|disease_id|disease_name|effect_text|disease_description

-- diseases.csv
CREATE TABLE diseases (
uniprot TEXT,
-- tmp_refs INTEGER,
disease_name TEXT,
effect_text TEXT,
disease_description TEXT,
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);


-- kinase_function_uniprotxml.csv
CREATE TABLE kin_function (
uniprot TEXT,
prot_function TEXT,
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);




-- kinase_reactions_uniprotxml.csv
CREATE TABLE reactions(
uniprot TEXT,
react_id INTEGER,
reaction_text TEXT,
PRIMARY KEY(uniprot, react_id),
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);



-- kinase_references.csv
-- CREATE TABLE references_full(
-- uniprot TEXT,
-- reference_id INTEGER,
-- pub_type TEXT,
-- pub_date TEXT,
-- pub_name TEXT,
-- pub_vol TEXT,
-- pub_pages_first TEXT,
-- pub_pages_last TEXT,
-- pub_title TEXT,
-- aut_text TEXT,
-- pubmedid TEXT,
-- doi TEXT,
-- PRIMARY KEY(uniprot, reference_id)
-- );



-- kin_subcell_loc_text.csv
CREATE TABLE subcell_location_text (
uniprot TEXT,
-- subcell_aditional_text_refs INTEGER,
subcell_aditional_text TEXT,
FOREIGN KEY(uniprot)  REFERENCES kinase_info(uniprot_id)
);


-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

--- de aqui para abajo @@@
-- kin_subcell_loc.csv
CREATE TABLE subcell_location (
uniprot TEXT,
subcell_location TEXT,
-- subcell_refs INTEGER,
FOREIGN KEY(uniprot)  REFERENCES kinase_info(uniprot_id)
);

-- Kinase_Substrate_Dataset.csv
CREATE TABLE kinase_substrate (
gene TEXT,
kinase TEXT,
kin_acc_id TEXT,
substrate TEXT,
sub_gene_id INTEGER,
sub_acc_id TEXT,
sub_gene TEXT,
sub_mod_rsd TEXT,
site_grp_id INTEGER,
site_7_aa TEXT,
prot_domain TEXT,
in_vivo_rxn TEXT,
in_vitro_rxn TEXT,
FOREIGN KEY (sub_acc_id) REFERENCES kinase_info(uniprot_id)
);



-- kinase_modified_residues_references_ref_numbers.csv
-- CREATE TABLE phosphosites_references (
--  uniprot_id TEXT,
--  residue_position INTEGER,
-- reference_id INTEGER,
--  --  FOREIGN KEY (uniprot_id, residue_position) REFERENCES phosphosites (uniprot_id, residue_position),
--  FOREIGN KEY(uniprot_id,reference_id ) REFERENCES references_full(uniprot, reference_id)
-- );


-- Group-project/csv_tables/kinases/isoforms/kinase_gral_info.csv
-- uniprot_id|full_prot_name|reverse|chromosome|start_gene_coord|genom_end_coord|sequence

CREATE TABLE isoforms_info (
uniprot_id TEXT,
full_prot_name TEXT,
reverse TEXT,
chromosome TEXT,
start_gene_coord INT,
genom_end_coord INT,
prot_sequence TEXT,
FOREIGN KEY(uniprot_id) REFERENCES  isoforms(isoforms)
);

-- csv_tables/ncbi_chromosomes
CREATE TABLE ncbi_chrom_id (
chr TEXT PRIMARY KEY,
ncbi_id TEXT
);


--------------------------------------------------
--------------------------------------------------
-- inhibitors

-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_table_unique_information.csv
-- inn_name,phase,mw,image_url,canonical_smiles,inchikey
CREATE TABLE inhibitors_gral_info (
inn_name TEXT PRIMARY KEY,
phase REAL,
mw REAL,
image_url TEXT ,
canonical_smiles TEXT,
inchikey TEXT
);

-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_synonims_table.csv
-- inn_name,synonyms
CREATE TABLE inhibitors_synonims (
inn_name  TEXT,
synonyms TEXT,
FOREIGN KEY  (inn_name) REFERENCES inhibitors_gral_info(inn_name)
);


-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitor_kinase_families.csv
-- inn_name,kinase_families
CREATE TABLE inhibitors_kin_family(
inn_name TEXT,
kinase_families TEXT,
FOREIGN KEY  (inn_name) REFERENCES inhibitors_gral_info(inn_name),
FOREIGN KEY(kinase_families) REFERENCES families(family_abbreviation)
);



-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/Inhibitors_pdbID.csv
-- inn_name,pdbid
CREATE TABLE inhibitors_pdbid (
inn_name TEXT ,
pdbid TEXT,
FOREIGN KEY  (inn_name) REFERENCES inhibitors_gral_info(inn_name)
);

-- /homes/dtg30/Desktop/group_proj/venv/src/Group-project/csv_tables/inhibitors/cleaned/inhibitors_target.csv
-- inn_name,targets
CREATE TABLE inhibitors_targets (
inn_name TEXT,
targets TEXT,
FOREIGN KEY  (inn_name) REFERENCES inhibitors_gral_info(inn_name)
);



-- Group-project/csv_tables/kinases/kinase_alternative_names.csv
CREATE TABLE kinase_alternative_names (
uniprot_id TEXT,
name TEXT,
short TEXT,
FOREIGN KEY(uniprot_id)  REFERENCES kinase_info(uniprot_id)
);


