from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect
import os
import queries



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
    # get arguments from url
    requested_name = request.values['search']
    requested_type = request.values['type']

    # modify the results according if the string was found as the name of the inhibitor
    #   or as a synonym
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

    # get required information from the database
    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], '*','inhibitors_gral_info',\
                                    ' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, basic.uniprot_id',\
                                  'inhibitors_targets inh LEFT JOIN  basic_info basic ON basic.gene = inh.targets',\
                                  ' inn_name LIKE "{}"'.format(inhib_name))
    synonyms = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    pdbid = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                      ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
    families = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])

    # pass the values to a dictionary, to be used inside the template
    context = {'gral_info':gral_info, 'targets': targets, 'synonyms':synonyms, 'pdbid':pdbid, 'families':families}

    return render_template('inhibitor_data.html', context=context)



@inhibitor_blueprint.route('/basic_info/<inhib_name>')
def inhibitor_data_information_iframe(inhib_name):
    ''' This is used to generate the <iframe> inside a modal on the kinase_data page, it
        gets and displays some basic information of the inhibitor
    '''

    # get information from the database
    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], '*','inhibitors_gral_info',' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, basic.uniprot_id',\
                                  'inhibitors_targets inh LEFT JOIN basic_info basic ON basic.prot_name = inh.targets',\
                                  ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'targets'])
    synonyms = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    families = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])

    context = {'gral_info':gral_info, 'targets': ', '.join(targets), 'synonyms': ', '.join(synonyms),  'families': ', '.join(families)}

    return render_template('inhibitor_data_information.html', context=context)



@inhibitor_blueprint.route('/<inhib_name>.json', methods=['GET'])
def inhibitor_json(inhib_name):
    ''' Returns a json page with the information requested on the url, useful as an api service.
        Only works for one inhibitor each request.
    '''

    # get column argument from the url, default value is ''
    requested_columns = request.args.get('columns','')
    requested_columns = requested_columns.rstrip().lower()

    # IF: the name of the inhibitor exists and is unique do the rest, else display null page
    if queries.query_is_unique(inhibitor_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info', \
                               'inn_name LIKE "{0}"'.format(inhib_name)):
        # get requested columns
        requested_columns_list = requested_columns.split(',')
        # if user is requesting all or didn't put anything, then get all avaliable  information
        if (('all' in requested_columns) or (requested_columns == '')):

            gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                            'inhibitors_gral_info',\
                                            ' inn_name LIKE "{}"'.format(inhib_name))
            targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets,  basic.uniprot_id',\
                                      'inhibitors_targets inh LEFT JOIN basic_info basic ON basic.prot_name = inh.targets',\
                                      ' inn_name LIKE "{}"'.format(inhib_name))
            synonyms_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'synonyms','inhibitors_synonims',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
            pdbid_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'pdbid','inhibitors_pdbid',\
                                              ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
            families_list = list(queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'kinase_families','inhibitors_kin_family',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])

            # generate dictionary containg all the information
            output_dict = {'general_info':gral_info.to_dict('index')[0],
                           'targets':targets.to_dict('index'),
                           'synonyms':synonyms_list, 'Protein_databank_id':pdbid_list,
                           'families':families_list
                            }
            # return the json
            return jsonify(output_dict)
        else:
            # if user requested particular information, get that from the database

            output_dict = {}
            for column  in requested_columns_list:
                # for loop to get required information. If requested, add it to the output dicionary
                if ('info' in column):
                    gral_info = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                                    'inhibitors_gral_info',\
                                                    ' inn_name LIKE "{}"'.format(inhib_name))
                    output_dict['general_info'] = gral_info.to_dict('index')[0]
                    
                elif ('targets' in column):
                    targets = queries.select_gral(inhibitor_blueprint.config['DATABASE'], 'inh.targets, basic.uniprot_id',\
                          'inhibitors_targets inh LEFT JOIN basic_info basic ON basic.prot_name = inh.targets',\
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

            # return json of the dictionary that contains the required information
            return jsonify(output_dict)
    else:
        # if inhibitor is not unique, then return ''
        return ''
  
