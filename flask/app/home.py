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

main_blueprint = Blueprint(
    'main', __name__,
    template_folder = 'templates/main',
)

main_blueprint.config = {}

@main_blueprint.record
def record_params(setup_state):
  app = setup_state.app
  main_blueprint.config = dict([(key,value) for (key,value) in app.config.items()])



#### home
@main_blueprint.route('/', methods=['GET','POST'])
def home():
    search_form = forms.Search_string()
    uploadfile_form = forms.UploadForm()
    phosphosite_form = forms.Phosphosite()
    inhibitor_form = forms.Inhibitors()
    context = {'uploadfile_form': uploadfile_form, 
               'search_form': search_form, 'phosphosite_form': phosphosite_form,
               'inhibitor_form':inhibitor_form}

    # kinase search form
    if request.method == 'POST' and search_form.search_string.data and search_form.validate_on_submit():
        requested_name = request.values['search_string']
        requested_name = requested_name.rstrip()
        if queries.query_is_unique(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info', \
                                   'uniprot_id LIKE "{0}" OR name_human LIKE "{0}" OR prot_name LIKE "{0}"'.format(requested_name)):
            uniprot_id = queries.select_gral(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info',\
                                             'uniprot_id LIKE "{0}" OR name_human LIKE "{0}" OR prot_name LIKE "{0}"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase.kinase_data',kin_name= uniprot_id))
        elif queries.query_is_unique(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info', \
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name)): # modify: get list of kinases
            uniprot_id = queries.select_gral(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info', \
                                             'uniprot_id LIKE "%{0}%"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase.kinase_data',kin_name= uniprot_id))
        elif queries.query_n_results(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info',\
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name))>1: # modify: get list of kinases
            return redirect(url_for('kinase.kinase_search_result', search=requested_name))
        else:
            text_flash = "'"+ str(requested_name) + "'" + " not in our database"
            flash(text_flash, 'kinase')
            return render_template('home.html', context = context)
        # modify:
        #  - if name is equal to a uniprot identifier redirect to /kinase/uniprotid
        # - else if requested_name correspond to one gene, redirect to the kinase/uniprotid
        # - else if requested_name one name nor to one gene, make a less restrictive querry and redirect to kinase_search_results


    # inhibitors form
    if request.method == 'POST' and inhibitor_form.data and inhibitor_form.validate_on_submit():
        inhib_requested = request.values['inhibitor_name']
        inhib_requested = inhib_requested.rstrip()
        if queries.query_is_unique(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info', \
                                   'inn_name LIKE "{0}"'.format(inhib_requested)):
            inn_name = queries.select_gral(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info',\
                                             'inn_name LIKE "{0}"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor.inhibitor_data',inhib_name= inn_name))
        elif queries.query_is_unique(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "{0}"'.format(inhib_requested)): # modify: get list of kinases
            inn_name = queries.select_gral(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "{0}"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor.inhibitor_data',inhib_name= inn_name))
        if queries.query_is_unique(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info', \
                                   'inn_name LIKE "%{0}%"'.format(inhib_requested)):
            inn_name = queries.select_gral(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_gral_info',\
                                             'inn_name LIKE "%{0}%"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor.inhibitor_data',inhib_name= inn_name))
        elif queries.query_is_unique(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested)): # modify: get list of kinases
            inn_name = queries.select_gral(main_blueprint.config['DATABASE'], 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor.inhibitor_data',inhib_name= inn_name))
        elif queries.query_n_results(main_blueprint.config['DATABASE'],  'inn_name','inhibitors_gral_info',\
                                     'inn_name LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('inhibitor.inhibitor_search_result', search=inhib_requested,type='inn'))
        elif queries.query_n_results(main_blueprint.config['DATABASE'],  'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('inhibitor.inhibitor_search_result', search=inhib_requested,type='syn'))

        else:
            text_flash = "'"+ str(inhib_requested) + "'" + " not in our database"
            flash(text_flash, 'inhibitor')
            return render_template('home.html', context = context)

    #phosphosite form
    if request.method == 'POST' and phosphosite_form.chromosome.data and phosphosite_form.validate_on_submit():
        chrom_requested = request.values['chromosome']
        chrom_requested = chrom_requested.rstrip()


    # uploaded file form
    if request.method == 'POST' and uploadfile_form.uploaded_file.data and uploadfile_form.validate_on_submit():
        file = request.files['uploaded_file']
        p_val_threshold = request.values['threshold_pval']
        threshold_foldchange = request.values['threshold_foldchange']
        threshold_cv_treatment = request.values['cv_treatment_threshold']
        random_name = str(random()).split('.')[1] #random number
        filename = random_name + '.tsv'
        session['tmp_upload_file'] = filename
        file.save(os.path.join(main_blueprint.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(main_blueprint.config['UPLOAD_FOLDER'], filename)

        return redirect(url_for('phosphoproteomics.phosphoproteomics',fc=threshold_foldchange,pv=p_val_threshold, cvt=threshold_cv_treatment), code = 307)

    
    return render_template('home.html', context = context)

