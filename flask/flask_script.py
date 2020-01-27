# lines to change signaled with "modify:   !!"

from flask import Flask, render_template, redirect, request, url_for, session
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



import pathlib
import re
path = str(pathlib.Path(__file__).parent.absolute())
DATABASE = re.sub(r'flask$','csv_tables/kinase_project.db',path)


#UPLOAD_FOLDER = '/home/daniel/Escritorio/uk/group_proj2/upload'
UPLOAD_FOLDER = '/homes/dtg30/Desktop/group_proj_2/'

app = Flask(__name__)
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER



#### home
@app.route('/', methods=['GET','POST'])
def home():
    search_form = forms.Search_string()
    uploadfile_form = forms.UploadForm()
    inhibitor_form = forms.Inhibitor_used()
    threshold_pval_form = forms.Threshold_pval()
    threshold_foldchange_form = forms.Threshold_foldchange()
    context = {'uploadfile_form': uploadfile_form, 'inhibitor_form': inhibitor_form,
               'threshold_foldchange_form': threshold_foldchange_form,
               'threshold_pval_form':threshold_pval_form,
               'search_form': search_form}

    if request.method == 'POST' and search_form.validate_on_submit():
        requested_name = request.values['search_string']
        requested_name = requested_name.rstrip()
        if queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info', 'uniprot_id LIKE "{}"'.format(requested_name)) : # modify: get list of kinases
            url = '/kinase/' + requested_name
            return redirect(url)
        elif queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info',  'name_human LIKE "{}"'.format(requested_name)) :
            uniprot_id = queries.select_gral(DATABASE, 'uniprot_id','kinase_info', 'name_human LIKE "{}"'.format(requested_name)).iloc[0,0]
            print(uniprot_id)
            url = '/kinase/' + uniprot_id
            return redirect(url)
        elif queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info', 'prot_name', requested_name) : # modify: get list of kinases
            uniprot_id = queries.select_gral(DATABASE, 'uniprot_id','kinase_info', 'prot_name LIKE "{}"'.format(requested_name)).iloc[0,0]
            print(uniprot_id)
            url = '/kinase/' + uniprot_id
            return redirect(url)
        if queries.query_n_results(DATABASE, 'uniprot_id','kinase_info', 'prot_name', requested_name) == 0: # modify: get list of kinases
            return render_template('home.html', context = context)




#        elif:

        # modify:
        #  - if name is equal to a uniprot identifier redirect to /kinase/uniprotid
        # - else if requested_name correspond to one gene, redirect to the kinase/uniprotid
        # - else if requested_name one name nor to one gene, make a less restrictive querry and redirect to kinase_search_results
        url = '/kinase/' + requested_name
        print(url)
        return redirect(url)

    if request.method == 'POST' and uploadfile_form.validate_on_submit():
        file = request.files['uploaded_file']
        p_val_threshold = request.values['threshold_pval']
        threshold_foldchange = request.values['threshold_foldchange']
        inhibitor = request.values['inhibitor']
        random_name = str(random()).split('.')[1] #random number
        filename = random_name + '.tsv'
        session['tmp_upload_file'] = filename
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        return redirect(url_for('phosphoproteomics',inh=inhibitor,fc=threshold_foldchange,pv=p_val_threshold, ), code = 307)

    return render_template('home.html', context = context)



#### kinases


@app.route('/kinase/results')
def kinase_search_result():
    # method: get , add filter
    return render_template('kinase_search_results.html')

@app.route('/kinase/<kin_name>')
def kinase_data(kin_name):
    # modify: method: post, add filter !!
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info', "uniprot_id LIKE '{}'".format(kin_name)) : # modify: get list of kinases
        kin_name = kin_name
        gral_info = queries.select_gral(DATABASE, '*', 'kinase_info', 'uniprot_id LIKE "{}"'.format(kin_name))
        isoforms_list = list(queries.select_gral(DATABASE, 'isoform', 'isoforms', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'isoform'])
        isoforms = ', '.join(isoforms_list)
        function_list = list(queries.select_gral(DATABASE, 'prot_function', 'kin_function', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'prot_function'])
        # modify: add layers of information
        reactions_list = list(queries.select_gral(DATABASE, 'reaction_text', 'reactions', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
        cell_loc_list= list(queries.select_gral(DATABASE, 'subcell_location', 'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])
        cell_loc_add_text_list = list(queries.select_gral(DATABASE, 'subcell_aditional_text', 'subcell_location_text', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_aditional_text'])
        diseases = queries.select_gral(DATABASE, 'DISTINCT disease_name, effect_text, disease_description', 'diseases', 'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
        context = {'kin_name':kin_name, 'gral_info': gral_info, 'isoforms': isoforms, 'function_list': function_list,
                   'reactions_list': reactions_list, 'cell_loc_list': cell_loc_list,
                   'cell_loc_add_text_list': cell_loc_add_text_list,
                   'diseases': diseases}
        return render_template('kinase_data.html', context = context)
    else:
        return 'not found'



#### inhibitors


@app.route('/inhibitor/search')
def inhibitor_search_result():
    # method: get, add filter
    return render_template('inhibitor_search_results.html')

@app.route('/inhibitor/<inhib_name>')
def inhibitor_data(inhib_name):
    # method: post, add filter
    return render_template('inhibitor_data.html')


### phosphoproteomics


@app.route('/phosphoproteomics', methods = ['GET','POST'])
def phosphoproteomics():
    context = {}
    if request.method == 'POST':
        tmp_file_name = session['tmp_upload_file'] # get name of the file
        session['tmp_upload_file'] = ''
        inhibitor = request.values['inh']
        fold_threshold = request.values['fc']
        pval_threshold = request.values['pv']
        tmp_file_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_file_name)
        try:
            ddf = phosphoproteomics_script.change_column_names(tmp_file_path, inhibitor)
            results_volcano = phosphoproteomics_script.volcano(ddf,pval_threshold, fold_threshold )
#            tims_function = phosphoproteomics_script.name(ddf, ...) # modify
        except:
            return 'Impossible to calculate, something wrong in the input values. <a href="/"> Go back </a>'
        context['volcano'] = results_volcano
        context['fold_threshold'] = fold_threshold
        context['pval_threshold'] = pval_threshold
    return render_template('phosphoproteomics.html', context = context)

### Documentation

@app.route('/documentation/general')
def documentation_general():
    return render_template('documentation_general.html')

@app.route('/documentation/api')
def documentation_api():
    return render_template('documentation_api.html')

@app.route('/documentation/stats')
def documentation_stats():
    x = [1,2,3]
    return render_template('documentation_stats.html', x = x, y = x)





### API

# return fasta file
# !! this does the api thing with .fasta return
# this could be done also to generate a csv

@app.route('/kinase/<kin_name>.fasta')
def fasta_protein(kin_name):
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info','uniprot_id  LIKE "{}"'.format(kin_name)):
        text = queries.select_gral(DATABASE, 'prot_sequence','kinase_info', 'uniprot_id  LIKE "{}"'.format(kin_name)) 
        sequence = text.loc[0,'prot_sequence']
        divide_each = 40  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + '| length: ' + str(seq_size) #modify: create header !!
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}
    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}



@app.route('/kinase/gene/<kin_name>.fasta')
def fasta_gene(kin_name):
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info','uniprot_id  LIKE "{}"'.format( kin_name)) : # modify: get list of kinases
        text = queries.select_gral(DATABASE, 'chromosome, genome_sequence, reverse, ensembl_gene_id, genome_starts, genome_ends','kinase_info', 'uniprot_id  LIKE "{}"'.format( kin_name))
        sequence = text.loc[0,'genome_sequence']
        divide_each = 40  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse'] + \
            ', Ensembl ID: ' + str( text.loc[0,'ensembl_gene_id']) + ', Start: ' + str(text.loc[0,'genome_starts']) + \
            ', Ends: ' + str(text.loc[0,'genome_ends'])
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}
    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}

# return
# this returns json !!
@app.route('/kinase/<kin_name>.json', methods=['GET'])
def api_all(kin_name):
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info','uniprot_id', kin_name) : # modify: get list of kinases
        q_output = queries.select_gral(DATABASE, '*','kinase_info', 'uniprot_id  LIKE "{}"'.format(kin_name))
        output_dict = [q_output.to_dict('index')[0]]
        return jsonify(output_dict)
    else:
        return ''


if __name__ == '__main__':
    app.run(debug=True, port = 8000)
