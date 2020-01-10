import retrieve_from_xml_uniprotapi as xmluniprot

uniprot = 'P31749'

data = xmluniprot.get_uniprot_data(uniprot)


grainfo = xmluniprot.get_gralinfo(data)

alt_prot = xmluniprot.get_alternative_prot_names(data)

ensembl_id = xmluniprot.get_ensembl_geneid(data)

modifres_gralinfo = xmluniprot.modres_gral_info(data)

modres_refs = xmluniprot.modres_references(data)
