from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
#from flask_csv import send_csv
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect
import os
import pandas as pd
import numpy as np
import phosphoproteomics_script
import sqlite3
from random import random
import queries
import divide_sequences
import add_pubmed_link


kinase_blueprint = Blueprint(
    'kinase', __name__,
    template_folder = 'templates/kinase',
    url_prefix = '/kinase'
)


@kinase_blueprint.record
def record_params(setup_state):
  app = setup_state.app
  kinase_blueprint.config = dict([(key,value) for (key,value) in app.config.items()])


#### kinases


@kinase_blueprint.route('/results')
def kinase_search_result():
    '''used when the user search for a kinase and more than 1 coincides with their search string'''
    # get arguments from url
    requested_name = request.values['search']

    # get results for the search
    search_results = queries.select_gral(kinase_blueprint.config['DATABASE'], 'uniprot_id, name_human, chromosome, fasd_name, ensembl_gene_id','kinase_info', 'uniprot_id LIKE "%{0}%"'.format(requested_name))
    
    context = {'search_results':search_results}

    return render_template('kinase_search_results.html', context = context)


@kinase_blueprint.route('/<kin_name>')
def kinase_data(kin_name):
    '''used when the user search for a kinase and only 1 coincides with their search string'''

    # check that the kinase is unique
    if queries.query_is_unique(kinase_blueprint.config['DATABASE'], 'uniprot_id', 'kinase_info',\
                               "uniprot_id LIKE '{}'".format(kin_name)) :
        # get information for that particular kinase stored in the database
        kin_name = kin_name
        gral_info = queries.select_gral(kinase_blueprint.config['DATABASE'], '*', 'kinase_info',\
                                        'uniprot_id LIKE "{}"'.format(kin_name))
        alternative_names = queries.select_gral(kinase_blueprint.config['DATABASE'], 'name, short',\
                                                'kinase_alternative_names' ,'uniprot_id LIKE "{}"'.format(kin_name))
        
        isoforms = list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'isoform', 'isoforms', \
                                            'uniprot LIKE "{}"'.format(kin_name)).loc[:,'isoform'])

         # inhibitors affecting that kinase
        inhibitors = list(queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                              'inn_name' ,'inhibitors_targets', \
                                              'targets LIKE "{}"'.format(gral_info.loc[0,'prot_name']) ).loc[:,'inn_name'])

        # kinase function, text needs to be modified to add hyperlinks to referenced pubmed 
        function_list_tmp = list(queries.select_gral(kinase_blueprint.config['DATABASE'], \
                                                     'prot_function', 'kin_function', \
                                                     'uniprot LIKE "{}"'.format(kin_name)).loc[:,'prot_function'])
        function_list =[]
        for function in function_list_tmp:
            function_list.append(add_pubmed_link.pubmed_link(function))

        reactions_list = list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'reaction_text', 'reactions',\
                                                  'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
        cell_loc_list= list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'DISTINCT subcell_location', \
                                                'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])

        targets = queries.select_gral(kinase_blueprint.config['DATABASE'], \
                                      'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa', 'kinase_substrate', \
                                      'kin_acc_id LIKE "{}"'.format(kin_name))

        phosphosites = queries.select_gral(kinase_blueprint.config['DATABASE'], \
                                           'residue_position, modif, type_modif, genom_begin, genom_end', \
                                           'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))

        # needs to be modified, since the text can contain pubmed references
        cell_loc_add_text_list_tmp = list(queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                                              'subcell_aditional_text', 'subcell_location_text',\
                                                              'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_aditional_text'])
        cell_loc_add_text_list = []
        for cell_loc in cell_loc_add_text_list_tmp:
            cell_loc_add_text_list.append(add_pubmed_link.pubmed_link(cell_loc))
            
        diseases = queries.select_gral(kinase_blueprint.config['DATABASE'], \
                                       'DISTINCT disease_name, effect_text, disease_description', 'diseases', \
                                       'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))

        # divide protein and gene sequences
        prot_seq_list = divide_sequences.divide_sequences(gral_info.loc[0,'prot_sequence'], 50,10)
        gene_seq_list = divide_sequences.divide_sequences(gral_info.loc[0,'genome_sequence'], 50, 10)


        context = {'kin_name':kin_name, 'gral_info': gral_info, 'isoforms': isoforms,
                   'alternative_names':alternative_names,'function_list': function_list,
                   'reactions_list': reactions_list, 'cell_loc_list': cell_loc_list,
                   'cell_loc_add_text_list': cell_loc_add_text_list,
                   'diseases': diseases,
                   'gene_seq_list':gene_seq_list, 'prot_seq_list': prot_seq_list,
                   'targets':targets, 'phosphosites':phosphosites,
                   'inhibitors': inhibitors
                   }
        return render_template('kinase_data.html', context = context)
    else:
        # if protein not unique or not found
        return 'not found'


########################
# genome viewer
#

@kinase_blueprint.route('genome_viewer/<uniprot_id>')
def genome_viewer(uniprot_id):
    '''used to generate a page that is used as an <iframe/> inside the kinase_data.html.
     Use to generate the url definition of the ncbi genome viewer, to be able to embed it later as an iframe
    '''
    # get phosphosites
    phosphosites = queries.select_gral(kinase_blueprint.config['DATABASE'], \
                                       'residue_position, modif, type_modif, genom_begin, genom_end', \
                                       'phosphosites', 'uniprot_id LIKE "{}"'.format(uniprot_id))
    # and 
    gral_info = queries.select_gral(kinase_blueprint.config['DATABASE'], 'kin.reverse, ncbi.ncbi_id',\
                                    'kinase_info kin LEFT JOIN ncbi_chrom_id ncbi ON ncbi.chr = kin.chromosome',\
                                    'kin.uniprot_id LIKE "{}"'.format(uniprot_id))
    chromosome_ncbi = gral_info.loc[0,'ncbi_id']

    genom_browser_markers_list = []
    color = '0040FF' #blue
    genom_browser_v = ''

    # check if the gene is in the reverse strand, and interchange start and end
    if (gral_info.loc[0,'reverse']  == 'True'):
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_end']) +':' + \
            str(phosphosites.loc[i,'genom_begin']) + '|' + str(phosphosites.loc[i,'residue_position']) +' '+\
            phosphosites.loc[i,'type_modif']+ '|' + color) 
        genom_browser_v = str(phosphosites['genom_end'].min()-20)  + ':' +  str(phosphosites['genom_begin'].max() +20)

    else:
        # IF not in the reverse (is in the forward), no need to interchange start and end
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_begin']) +':' + \
            str(phosphosites.loc[i,'genom_end']) + '|' + str(phosphosites.loc[i,'residue_position']) + '|' + color)
        genom_browser_v = str(phosphosites['genom_begin'].min()-20) + ':' + str(phosphosites['genom_end'].max() +20)

    # string of markers
    genom_browser_markers = ','.join(genom_browser_markers_list)

    # generate url for ncbi genome viewer
    genom_browser_src ="?embedded=true&id="+chromosome_ncbi+"&tracks=[key:sequence_track,name:T16507,display_name:Sequence,id:T16507,dbname:SADB,annots:NA000001672.2,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:six_frames_translation,name:T11044,display_name:Six-frame translations,id:T11044,dbname:GenBank,annots:Six-frame translation,ShowOption:All,OrfThreshold:20,HighlightCodons:true,AltStart:false,shown:true,order:21][key:gene_model_track,name:T2000262,display_name:Genes\, NCBI Homo sapiens Annotation Release 105.20190906,id:T2000262,dbname:SADB,annots:NA000229419.1,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:22]&assm_context=GCF_000001405.13&mk="+genom_browser_markers+"&v="+genom_browser_v+"&c=ffff99&select=gi|224589800-0017f72a-001844c0-010a-ff7f5665-NA000229419.1;&slim=0&appname=pkinases"


    context = {'genom_browser_src':genom_browser_src}
    return render_template('genome_viewer.html', context=context)






### API

@kinase_blueprint.route('/<kin_name>.fasta')
def fasta_protein(kin_name):
    '''displays protein sequence in fasta format'''
    # check if requested kinase exists in the database and is unique
    if queries.query_is_unique(kinase_blueprint.config['DATABASE'], 'uniprot_id', 'kinase_info','uniprot_id  LIKE "{}"'.format(kin_name)):
        # get information from database
        text = queries.select_gral(kinase_blueprint.config['DATABASE'], 'prot_sequence, chromosome, reverse','kinase_info', 'uniprot_id  LIKE "{}"'.format(kin_name))
        sequence = text.loc[0,'prot_sequence']
        divide_each = 40  # size of each row
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            # generate header            
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse']
            text_out = '\n'.join([header, seq_out])

        # return page with fasta sequence
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}

    # IF: it is an isoform, do the same: display the protein sequence
    elif queries.query_is_unique(kinase_blueprint.config['DATABASE'], 'uniprot_id', 'isoforms_info','uniprot_id  LIKE "{}"'.format(kin_name)) :
        # get information of sequence and for the header
        text = queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                   'prot_sequence, chromosome, reverse','isoforms_info', \
                                   'uniprot_id  LIKE "{}"'.format( kin_name))
        sequence = text.loc[0,'prot_sequence']
        divide_each = 40  # characters per row
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            # generate header            
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse']
            text_out = '\n'.join([header, seq_out])
        # return fasta file
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}

    # if it was not a kinase or an isoform, or not found, return ''
    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}



@kinase_blueprint.route('/gene/<kin_name>.fasta')
def fasta_gene(kin_name):
    '''displays gene sequence in fasta format'''
    # check if requested kinase exists in the database and is unique
    if queries.query_is_unique(kinase_blueprint.config['DATABASE'], 'uniprot_id', 'kinase_info','uniprot_id  LIKE "{}"'.format( kin_name)) :
        # get information from database        
        text = queries.select_gral(kinase_blueprint.config['DATABASE'], 'chromosome, genome_sequence, reverse, ensembl_gene_id, genome_starts, genome_ends','kinase_info', 'uniprot_id  LIKE "{}"'.format( kin_name))
        sequence = text.loc[0,'genome_sequence']
        divide_each = 40  # size of each row
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            # generate header
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse'] + \
            ', Ensembl ID: ' + str( text.loc[0,'ensembl_gene_id']) + ', Start: ' + str(text.loc[0,'genome_starts']) + \
            ', Ends: ' + str(text.loc[0,'genome_ends'])
            text_out = '\n'.join([header, seq_out])
        # return fasta file            
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}
    else:
        # if it was not a kinase, return ''
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}


# this returns json !!
@kinase_blueprint.route('/<kin_name>.json', methods=['GET'])
def kinase_json(kin_name):
    ''' Returns a json page with the information requested on the url, useful as an api service.
        Only works for one inhibitor each request.
    '''
    # get column argument from the url, default value is ''
    requested_columns = request.args.get('columns','')
    requested_columns = requested_columns.rstrip().lower()

    # IF: the name of the kinase exists and is unique do the rest, else display null page
    if queries.query_is_unique(kinase_blueprint.config['DATABASE'], 'uniprot_id','kinase_info',\
                                   'uniprot_id LIKE "{0}"'.format(kin_name)):
        # get requested columns
        requested_columns_list = requested_columns.split(',')
        # IF: user is requesting all or didn't put anything, then get all avaliable  information
        if (('all' in requested_columns) or (requested_columns == '')):

            gral_info = queries.select_gral(kinase_blueprint.config['DATABASE'], '*','kinase_info', \
                                            'uniprot_id LIKE "{0}" '.format(kin_name))
            isoforms = queries.select_gral(kinase_blueprint.config['DATABASE'], 'kin.isoform, iso.prot_sequence',\
                                           'isoforms kin LEFT JOIN isoforms_info iso ON kin.isoform = iso.uniprot_id', \
                                           'kin.uniprot LIKE "{}"'.format(kin_name)) 
            inhibitors = queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                                  'inh.inn_name, inhinfo.phase, inhinfo.mw, inhinfo.canonical_smiles, inhinfo.inchikey' ,\
                                                  'kinase_info kin  INNER JOIN inhibitors_targets inh ON kin.prot_name = inh.targets INNER JOIN inhibitors_gral_info inhinfo ON inh.inn_name = inhinfo.inn_name',\
                                                  'kin.uniprot_id LIKE "{}"'.format(kin_name))
            reactions_list = list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'reaction_text', 'reactions',\
                                                      'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
            cell_loc_list= list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'DISTINCT subcell_location',\
                                                    'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])
            diseases = queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                           'DISTINCT disease_name, effect_text, disease_description', 'diseases',\
                                           'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
            targets = queries.select_gral(kinase_blueprint.config['DATABASE'], 'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa',\
                                          'kinase_substrate', 'kin_acc_id LIKE "{}"'.format(kin_name))
            phosphosites = queries.select_gral(kinase_blueprint.config['DATABASE'], 'residue_position, modif, type_modif, genom_begin, genom_end', \
                                               'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))
            # generate dictionary containg all the information
            output_dict = {'general_info':gral_info.to_dict('index')[0],
                           'isoforms': isoforms.to_dict('index'), 'inhibitors': inhibitors.to_dict('index'),
                           'reactions':reactions_list, 'cell_loc': cell_loc_list,
                           'diseases':diseases.to_dict('index'), 'targets':targets.to_dict('index'),
                           'phosphosites': phosphosites.to_dict('index')}
        else:
            # if requested particular information:
            output_dict = {}
            for column  in requested_columns_list:
                # for loop to get required information. If requested, add it to the output dicionary
                if ('info' in column):
                    gral_info = queries.select_gral(kinase_blueprint.config['DATABASE'], '*','kinase_info', \
                                                    'uniprot_id LIKE "{0}" '.format(kin_name))
                    output_dict['general_info'] = gral_info.to_dict('index')[0]

                elif ('isoforms' in column):
                    isoforms = queries.select_gral(kinase_blueprint.config['DATABASE'], 'kin.isoform, iso.prot_sequence',\
                                                   'isoforms kin LEFT JOIN isoforms_info iso ON kin.isoform = iso.uniprot_id', \
                                                   'kin.uniprot LIKE "{}"'.format(kin_name)) 
                    output_dict['isoforms'] = isoforms.to_dict('index')

                elif ('inhibitors' in column):
                    inhibitors = queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                                     'inh.inn_name, inhinfo.phase, inhinfo.mw, inhinfo.canonical_smiles, inhinfo.inchikey' ,\
                                                     'kinase_info kin  INNER JOIN inhibitors_targets inh ON kin.prot_name = inh.targets INNER JOIN inhibitors_gral_info inhinfo ON inh.inn_name = inhinfo.inn_name',\
                                                     'kin.uniprot_id LIKE "{}"'.format(kin_name))
                    output_dict['inhibitors'] = inhibitors.to_dict('index')

                elif ('reactions' in column):
                    reactions_list = list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'reaction_text', 'reactions',\
                                                              'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
                    output_dict['reactions'] = reactions_list

                elif ('cell_loc' in column):
                    cell_loc_list = list(queries.select_gral(kinase_blueprint.config['DATABASE'], 'DISTINCT subcell_location',\
                                                             'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])
                    output_dict['cell_loc'] = cell_loc_list

                elif ('diseases' in column):
                    diseases = queries.select_gral(kinase_blueprint.config['DATABASE'],\
                                                   'DISTINCT disease_name, effect_text, disease_description', 'diseases',\
                                                   'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
                    output_dict['diseases'] = diseases.to_dict('index')

                elif ('phosphosites_self' in column):
                    phosphosites = queries.select_gral(kinase_blueprint.config['DATABASE'], 'residue_position, modif, type_modif, genom_begin, genom_end', \
                                                       'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))
                    output_dict['phosphosites_self'] = phosphosites.to_dict('index')

                elif ('phosphosites_targets' in column):
                    targets = queries.select_gral(kinase_blueprint.config['DATABASE'], 'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa',\
                                                  'kinase_substrate', 'kin_acc_id LIKE "{}"'.format(kin_name))
                    output_dict['phosphosites_targets'] = targets.to_dict('index')
        # return json of the dictionary that contains the required information
        return jsonify(output_dict)
    

    else:
        # if kinase is not unique, then return ''
        return ''

