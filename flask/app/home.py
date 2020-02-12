from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
from flask import Response # for api: fasta , csv and so on
from flask_wtf import CsrfProtect
import os
from random import random
import queries

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

    # KINASE SEARCH FORM
    if request.method == 'POST' and search_form.search_string.data and search_form.validate_on_submit():

        # get values from the form
        requested_name = request.values['search_string']
        requested_name = requested_name.rstrip()

        # IF: requested value is unique redirect to info for that particular kinase (kinase_data.html)
        if queries.query_is_unique(main_blueprint.config['DATABASE'], 'kin.uniprot_id',\
                                   'kinase_info kin LEFT JOIN basic_info basic ON basic.uniprot_id = kin.uniprot_id',\
                                   'kin.uniprot_id LIKE "{0}" OR kin.name_human LIKE "{0}" OR basic.gene LIKE "{0}"'.format(requested_name)):
            uniprot_id = queries.select_gral(main_blueprint.config['DATABASE'], 'kin.uniprot_id',\
                                  'kinase_info kin LEFT JOIN basic_info basic ON basic.uniprot_id = kin.uniprot_id',\
                                  'kin.uniprot_id LIKE "{0}" OR kin.name_human LIKE "{0}" OR basic.gene LIKE "{0}"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase.kinase_data',kin_name= uniprot_id))
        elif queries.query_is_unique(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info', \
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name)): # modify: get list of kinases
            uniprot_id = queries.select_gral(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info', \
                                             'uniprot_id LIKE "%{0}%"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase.kinase_data',kin_name= uniprot_id))

        # IF: requested value is contained in several rows, redirect to result table (kinase_search_result.html)
        elif queries.query_n_results(main_blueprint.config['DATABASE'], 'uniprot_id','kinase_info',\
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name))>1: # modify: get list of kinases
            return redirect(url_for('kinase.kinase_search_result', search=requested_name))

        # IF: nothing was found, stay in the home, and tell the user that search string was not found
        #        in our database
        else:
            text_flash = "'"+ str(requested_name) + "'" + " not in our database"
            flash(text_flash, 'kinase')
            return render_template('home.html', context = context)

    # INHIBITORS FORM
    # validate on submit, assure that contains data
    if request.method == 'POST' and inhibitor_form.inhibitor_name.data and inhibitor_form.validate_on_submit():

        # get values from FORM
        inhib_requested = request.values['inhibitor_name']
        inhib_requested = inhib_requested.rstrip()

        # IF: requested value is unique redirect to info for that particular inhibitor (inhibitor_data.html)
        #    check in inhibitor name and synonyms, first with a strict match the with soft matching
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
        
        # IF: requested value is contained in several rows, redirect to result table (inhibitor_search_result.html)
        elif queries.query_n_results(main_blueprint.config['DATABASE'],  'inn_name','inhibitors_gral_info',\
                                     'inn_name LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('inhibitor.inhibitor_search_result', search=inhib_requested,type='inn'))
        elif queries.query_n_results(main_blueprint.config['DATABASE'],  'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('inhibitor.inhibitor_search_result', search=inhib_requested,type='syn'))

        # IF: nothing was found, stay in the home, and tell the user that search string was not found
        #        in our database
        else:
            text_flash = "'"+ str(inhib_requested) + "'" + " not in our database"
            flash(text_flash, 'inhibitor')
            return render_template('home.html', context = context)

    # PHOSPHOSITE FORM:
    if request.method == 'POST' and phosphosite_form.chromosome.data and phosphosite_form.validate_on_submit():
        # extract values from the FORM
        chrom_requested = request.values['chromosome']
        genome_start = int(request.values['genomic_loc_start'])
        genome_ends = int(request.values['genomic_loc_end'])

        if (genome_ends <= genome_start):
          flash('start should be < than end point','genomic_location')
          return render_template('home.html', context = context)

        # get chromosome number and ncbi chromosome id
        chrom_number,crhom_ncbi = chrom_requested.split(',')

        return redirect(url_for('genomic_browser.genome_viewer_chrom',chromosome=chrom_number,ncbi=crhom_ncbi,gnm_start=genome_start, gnm_end=genome_ends), code = 307)


    # UPLOADED FILE FORM:
    # validate and check that there is an uploaded file
    if request.method == 'POST' and uploadfile_form.uploaded_file.data and uploadfile_form.validate_on_submit():

        # extract values from the FORM
        file = request.files['uploaded_file']
        p_val_threshold = request.values['threshold_pval']
        threshold_foldchange = request.values['threshold_foldchange']
        threshold_cv_treatment = request.values['cv_treatment_threshold']

        # save file with a random name, to prevent conflicts
        random_name = str(random()).split('.')[1] #random number
        filename = random_name + '.tsv'
        session['tmp_upload_file'] = filename
        file.save(os.path.join(main_blueprint.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(main_blueprint.config['UPLOAD_FOLDER'], filename)

        return redirect(url_for('phosphoproteomics.phosphoproteomics',fc=threshold_foldchange,pv=p_val_threshold, cvt=threshold_cv_treatment), code = 307)

    
    return render_template('home.html', context = context)

