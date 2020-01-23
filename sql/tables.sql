
-- kinase_general_information.csv
CREATE TABLE kinase_info (
uniprot_id	VARCHAR(50) PRIMARY KEY,
full_prot_name	VARCHAR(100),
reverse	VARCHAR(10),
chromosome	VARCHAR(2),
prot_sequence	TEXT,
family	VARCHAR(10),
prot_name	VARCHAR(100),
name_human	VARCHAR(20),
mass	INT,
ensembl_gene_id	VARCHAR(20),
genome_starts	INT,
genome_ends	INT,
genome_sequence TEXT
FOREIGN KEY(family) REFERENCES families(family_abbreviation)
);

-- kinase_isoforms_uniprotxml.csv
CREATE TABLE isoforms (
uniprot VARCHAR(20),
isoform VARCHAR(20) PRIMARY KEY,
FOREIGN KEY(uniprot)  REFERENCES kinase_info(uniprot_id)
);

-- families.csv
CREATE TABLE families (
family_abbreviation VARCHAR(20) PRIMARY KEY,
residue VARCHAR(20),
family_name VARCHAR(100)
);

-- phosphosite_table.csv
CREATE TABLE phosphosites (
uniprot_id	VARCHAR(20),
residue_position	INT,
modif	VARCHAR(30),
type_modif	VARCHAR(50),
genom_begin	INT,
genom_end  INT
PRIMARY KEY(uniprot_id, residue_position),
FOREIGN KEY(uniprot_id) REFERENCES kinase_info(uniprot_id)
);

-- diseases.csv
CREATE TABLE diseases (
id INT PRIMARY KEY AUTOINCREMENT,
uniprot	VARCHAR(20),
tmp_refs	INT,
disease_name	VARCHAR(20),
effect_text	VARCHAR(400),
disease_description TEXT,
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);

-- kinase_function_uniprotxml.csv
CREATE TABLE kin_function (
id INT PRIMARY KEY AUTOINCREMENT,
uniprot	VARCHAR(20),
prot_function TEXT,
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);

-- kinase_function_refs_uniprotxml.csv
CREATE TABLE function_references (
id INT PRIMARY KEY AUTOINCREMENT,
uniprot	VARCHAR(20),
item INT,
FOREIGN KEY(uniprot) REFERENCES kin_function(uniprot)
);


-- kinase_reactions_uniprotxml.csv
CREATE TABLE reactions(
uniprot VARCHAR(20),
react_id	INT,
reaction_text VARCHAR(200),
PRIMARY KEY(uniprot, react_ref),
FOREIGN KEY(uniprot) REFERENCES kinase_info(uniprot_id)
);

-- kinase_reactions_refs_uniprotxml.csv
CREATE TABLE reactions_refences (
ref_item	INT,
uniprot	VARCHAR(20),
react_id INT,
PRIMARY KEY(uniprot,ref_item),
FOREIGN KEY (uniprot, react_id) REFERENCES reactions(uniprot, react_id)
);

-- kin_subcell_loc_text.csv
CREATE TABLE subcell_location_text (
id INT PRIMARY KEY AUTOINCREMENT,
uniprot	VARCHAR(20),
subcell_aditional_text_refs	INT,
subcell_aditional_text TEXT
FOREIGN KEY(uniprot)  REFERENCES kinase_info(uniprot_id),
FOREIGN KEY(uniprot, subcell_aditional_text_refs) REFERENCES references_full(uniprot, reference_id)
);

-- kin_subcell_loc.csv
CREATE TABLE subcell_location (
id INT PRIMARY KEY AUTOINCREMENT,
uniprot	VARCHAR(20),
subcell_location	VARCHAR(20),
subcell_refs INT
FOREIGN KEY(uniprot, subcell_refs) REFERENCES references_full(uniprot, reference_id)
);

-- Kinase_Substrate_Dataset.csv
CREATE TABLE kinase_substrate (
id INT PRIMARY KEY AUTOINCREMENT,
gene	VARCHAR(20),
kinase	VARCHAR(20),
kin_acc_id	VARCHAR(20),
substrate	VARCHAR(20),
sub_gene_id	INT,
sub_acc_id	VARCHAR(20),
sub_gene	VARCHAR(20),
sub_mod_rsd	VARCHAR(20),
site_grp_id	INT,
site_7_aa	VARCHAR(20),
prot_domain	VARCHAR(20),
in_vivo_rxn	VARCHAR(2),
in_vitro_rxn VARCHAR(20)
FOREIGN KEY (sub_acc_id) REFERENCES kinase_info(uniprot_id)
);


-- kinase_references.csv
CREATE TABLE references_full (
uniprot	VARCHAR(20),
reference_id	INT,
pub_type	VARCHAR(20),
pub_date	VARCHAR(20),
pub_name	VARCHAR(80),
pub_vol	VARCHAR(20),
pub_pages_first	VARCHAR(20),
pub_pages_last	VARCHAR(20),
pub_title	VARCHAR(300),
aut_text	TEXT,
pubmedid	VARCHAR(20),
doi VARCHAR(20)
PRIMARY KEY(uniprot, reference_id)
);
