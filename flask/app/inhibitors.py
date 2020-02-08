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



inhibitor_blueprint = Blueprint(
    'inhibitor', __name__,
    template_folder = 'templates/inhibitor',
    url_prefix = '/inhibitor'
)

@inhibitor_blueprint.record
def record_params(setup_state):
  app = setup_state.app
  inhibitor_blueprint.config = dict([(key,value) for (key,value) in app.config.items()])













#### inhibitors


@inhibitor_blueprint.route('/search')
def inhibitor_search_result():
    # method: get, add filter
    requested_name = request.values['search']
    requested_type = request.values['type']
    if (requested_type == 'inn'):
        search_results = queries.select_gral(inhibitor_blueprint.config['DATABASE'],  'inn_name, phase, mw','inhibitors_gral_info',\
                                     'inn_name LIKE "%{0}%"'.format(requested_name))
    elif (requested_type == 'syn'):
        search_results = queries.select_gral(inhibitor_blueprint.config['DATABASE'],  'inn_name, phase, mw','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(requested_name))
    context = {'search_results': search_results}
    return  render_template('inhibitor_search_results.html', context = context)

@inhibitor_blueprint.route('/<inhib_name>')
def inhibitor_data(inhib_name):
    # method: post, add filter
    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], '*','inhibitors_gral_info',\
                                    ' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, kin.uniprot_id',\
                                  'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                                  ' inn_name LIKE "{}"'.format(inhib_name))


    synonyms = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    pdbid = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                      ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
    families = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
    context = {'gral_info':gral_info, 'targets': targets, 'synonyms':synonyms, 'pdbid':pdbid, 'families':families}
    return render_template('inhibitor_data.html', context=context)



@inhibitor_blueprint.route('/basic_info/<inhib_name>')
def inhibitor_data_information_iframe(inhib_name):
    # method: post, add filter
    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], '*','inhibitors_gral_info',' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, kin.uniprot_id','inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',' inn_name LIKE "{}"'.format(inhib_name))


    synonyms = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    pdbid = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                      ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
    families = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
    context = {'gral_info':gral_info, 'targets': targets, 'synonyms':synonyms, 'pdbid':pdbid, 'families':families}
    return render_template('inhibitor_data_information.html', context=context)



@inhibitor_blueprint.route('/<inhib_name>.json', methods=['GET'])
def inhibitor_json(inhib_name):
    requested_columns = request.args.get('columns','')
    requested_columns = requested_columns.rstrip().lower()
    if queries.query_is_unique(inhibitor_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info', \
                               'inn_name LIKE "{0}"'.format(inhib_name)):
        requested_columns_list = requested_columns.split(',')
        if (('all' in requested_columns) or (requested_columns == '')):

            gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                            'inhibitors_gral_info',\
                                            ' inn_name LIKE "{}"'.format(inhib_name))
            targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, kin.uniprot_id',\
                                      'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                                      ' inn_name LIKE "{}"'.format(inhib_name))
            synonyms_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
            pdbid_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                              ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
            families_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
        
            output_dict = {'general_info':gral_info.to_dict('index')[0],
                           'targets':targets.to_dict('index'),
                           'synonyms':synonyms_list, 'Protein_databank_id':pdbid_list,
                           'families':families_list
                            }
            return jsonify(output_dict)
        else:
            output_dict = {}
            for column  in requested_columns_list:
                if ('info' in column):
                    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                                    'inhibitors_gral_info',\
                                                    ' inn_name LIKE "{}"'.format(inhib_name))
                    output_dict['general_info'] = gral_info.to_dict('index')[0]
                elif ('targets' in column):
                    targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, kin.uniprot_id',\
                          'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                          ' inn_name LIKE "{}"'.format(inhib_name))
                    output_dict['targets'] = targets.to_dict('index')
                elif ('synonyms' in column):
                    synonyms_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                     ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
                    output_dict['synonyms'] = synonyms_list
                elif ('pbid' in column):
                    pdbid_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                  ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
                    output_dict['Protein_databank_id'] = pdbid_list
                elif ('families' in column):
                    families_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                     ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
                    output_dict['families'] = families_list
            return jsonify(output_dict)
    else:
        return ''
  
