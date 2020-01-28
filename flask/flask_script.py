# lines to change signaled with "modify:   !!"

from flask import Flask, render_template, redirect, request, url_for, session, flash
import forms
#from flask_csv import send_csv
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect
from werkzeug.utils import secure_filename
import os
import pandas as pd
import numpy as np
import phosphoproteomics_script

#UPLOAD_FOLDER = '/home/daniel/Escritorio/uk/group_proj2/upload'
UPLOAD_FOLDER = '/homes/ta317/Desktop'


app = Flask(__name__)
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#

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
        requested_uniprot = request.values['search_string']
        # modify: if clause
        url = '/kinase/' + requested_uniprot
        print(url)
        return redirect(url)
    if request.method == 'POST' and uploadfile_form.validate_on_submit():
        file = request.files['uploaded_file']
        p_val_threshold = request.values['threshold_pval']
        threshold_foldchange = request.values['threshold_foldchange']
        inhibitor = request.values['inhibitor']
        print(inhibitor)
        #print(secure_filename(file.filename))
        filename = 'test.tsv' #!!! change this line
        session['tmp_upload_file'] = filename
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        return redirect(url_for('phosphoproteomics',inh=inhibitor,fc=threshold_foldchange,pv=p_val_threshold, ), code = 307)
    else:
        flash('You were successfully logged in') # ~
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
        isoforms = list(queries.select_gral(DATABASE, 'isoform', 'isoforms', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'isoform'])
#        isoforms = ', '.join(isoforms_list)

        function_list_tmp = list(queries.select_gral(DATABASE, 'prot_function', 'kin_function', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'prot_function'])
        function_list =[]
        for function in function_list_tmp:
            function_list.append(add_pubmed_link.pubmed_link(function))

        # modify: add layers of information
        reactions_list = list(queries.select_gral(DATABASE, 'reaction_text', 'reactions', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
        cell_loc_list= list(queries.select_gral(DATABASE, 'subcell_location', 'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])

        targets = queries.select_gral(DATABASE, 'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa', 'kinase_substrate', 'kin_acc_id LIKE "{}"'.format(kin_name))

        phosphosites = queries.select_gral(DATABASE, 'residue_position, modif, type_modif, genom_begin, genom_end', 'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))

        cell_loc_add_text_list_tmp = list(queries.select_gral(DATABASE, 'subcell_aditional_text', 'subcell_location_text', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_aditional_text'])
        cell_loc_add_text_list = []
        for cell_loc in cell_loc_add_text_list_tmp:
            cell_loc_add_text_list.append(add_pubmed_link.pubmed_link(cell_loc))

        diseases = queries.select_gral(DATABASE, 'DISTINCT disease_name, effect_text, disease_description', 'diseases', 'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
        prot_seq_list = divide_sequences.divide_sequences(gral_info.loc[0,'prot_sequence'], 50,10)
        gene_seq_list = divide_sequences.divide_sequences(gral_info.loc[0,'genome_sequence'], 50, 10)

        context = {'kin_name':kin_name, 'gral_info': gral_info, 'isoforms': isoforms,
                   'function_list': function_list,
                   'reactions_list': reactions_list, 'cell_loc_list': cell_loc_list,
                   'cell_loc_add_text_list': cell_loc_add_text_list,
                   'diseases': diseases,
                   'gene_seq_list':gene_seq_list, 'prot_seq_list': prot_seq_list,
                   'targets':targets, 'phosphosites':phosphosites}

    try:
        return render_template('kinase_data.html', context = context)
    except:
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
def phosphoproteomics():  #~~~
    context = {}
    if request.method == 'POST':
        tmp_file_name = session['tmp_upload_file'] # get something from request
        session['tmp_upload_file'] = ''
        inhibitor = request.values['inh']
        fold_threshold = request.values['fc']
        pval_threshold = request.values['pv']
        tmp_file_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_file_name)
        #modify : Query to get kinase substrate
        try:
            ddf = phosphoproteomics_script.change_column_names(tmp_file_path, inhibitor)
            results_volcano = phosphoproteomics_script.volcano(ddf,pval_threshold, fold_threshold )
            Z_score = phosphoproteomics_script.KSEA(ddf, kinase_substrate)
        except:
            return 'Impossible to calculate, something wrong in the input values. <a href="/"> Go back </a>'
        context['volcano'] = results_volcano
        context['fold_threshold'] = fold_threshold
        context['pval_threshold'] = pval_threshold
        context['non_identified'] = Z_score['non_identified']
        context['score'] = Z_score['score']
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
        text = queries.select_gral(DATABASE, 'prot_sequence, chromosome, reverse','kinase_info', 'uniprot_id  LIKE "{}"'.format(kin_name))
        sequence = text.loc[0,'prot_sequence']
        divide_each = 40  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse']
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}

    elif queries.query_is_unique(DATABASE, 'uniprot_id', 'isoforms_info','uniprot_id  LIKE "{}"'.format(kin_name)) : # modify: get list of kinases
        text = queries.select_gral(DATABASE, 'prot_sequence, chromosome, reverse','isoforms_info', 'uniprot_id  LIKE "{}"'.format( kin_name))
        sequence = text.loc[0,'prot_sequence']
        divide_each = 40  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse']
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}
    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}



@app.route('/kinase/gene/<kin_name>.fasta')
def fasta_gene(kin_name):
    if kin_name == '1': # modify: get list of kinases
        sequence = 'aoisjdoaisjdaoisdj' #modify: retrieve from database !! needs an if/elseto handle non existent
        divide_each = 10  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + '|' + str(seq_size) #modify: create header !!
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}

    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}


# return
# this returns json !!
@app.route('/kinase/<kin_name>.json', methods=['GET'])
def api_all(kin_name):
    # modify: add a if/else handler to check if that prot exist in database
    # modify: data should be retrieved from database. The following is just an example !!
    # modify: output_dict is a list of dictionaries
    arguments = request.args
    # arguments are obtained as tupples
    print(arguments)  # modify: get arguments, needs a handler, check minimum existent, filter output by arguments, ...
    output_dict =     [{'id': 2,
     'uniprot_accession': kin_name,
     'info1': 'text1',
     'info2': 'text2',
     'info3': 'text3'}]
    return jsonify(output_dict)


if __name__ == '__main__':
    app.run(debug=True, port = 8000)
