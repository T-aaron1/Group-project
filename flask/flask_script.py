# lines to change signaled with "modify:   !!"

from flask import Flask, render_template, redirect, request, url_for, session
import forms
#from flask_csv import send_csv
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect
from werkzeug.utils import secure_filename
import os
import pandas as pd
import numpy as np


#UPLOAD_FOLDER = '/home/daniel/Escritorio/uk/group_proj2/upload'
UPLOAD_FOLDER = '/homes/dtg30/Desktop/group_proj_2/'

app = Flask(__name__)
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#

#### home
@app.route('/', methods=['GET','POST'])
def home():
    uploadfile_form = forms.UploadForm()
    context = {'uploadfile_form': uploadfile_form}
    if request.method == 'POST' and uploadfile_form.validate_on_submit():
        file = request.files['uploaded_file']
        #print(secure_filename(file.filename))
        filename = 'test.tsv' #!!! change this line
        session['tmp_upload_file'] = filename
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        return redirect(url_for('phosphoproteomics'), code = 307)
    return render_template('home.html', context = context)



#### kinases

@app.route('/kinase')
def kinase_general():
    return render_template('kinase_general.html')

@app.route('/kinase/results')
def kinase_search_result():
    # method: get , add filter
    return render_template('kinase_search_results.html')

@app.route('/kinase/<kin_name>')
def kinase_data(kin_name):
    # modify: method: post, add filter !!
    kin_name = kin_name
    # modify: add layers of information
    sql_retrieve = 'this is sql result'
    context = {'kin_name':kin_name, 'sql': sql_retrieve}
    try:
        return render_template('kinase_data.html', context = context)
    except:
        return 'not found'



#### inhibitors
@app.route('/inhibitor')
def inhibitor_general():
    return render_template('inhibitor_general.html')

@app.route('/inhibitor/search')
def inhibitor_search_result():
    # method: get, add filter
    return render_template('inhibitor_search_results.html')

@app.route('/inhibitor/<inhib_name>')
def inhibitor_data(inhib_name):
    # method: post, add filter
    return render_template('inhibitor_data.html')


### phosphoproteomics


def volcano(file_path): #!! change this, mod threshold (from request), put in dif file, get min, and max for axis
    df = pd.DataFrame(pd.read_csv(file_path, sep='\t'))
    df.dropna(axis= 1, how='all', inplace = True)
    list_nulls = list(df[df.AZ20_fold_change.isnull()].Substrate)
    list_control_zero = list(df[df.control_mean == 0].Substrate)
    df = df[(df.AZ20_fold_change.notnull()) & (df.control_mean != 0) & (df.AZ20_fold_change != 0)]
    df['log_foldchange'] = np.log2(df['AZ20_fold_change'])
    df['log_pval'] = -np.log10(df['AZ20_p-value'])

    df =df[~(np.isinf(np.abs(df.log_foldchange)))& ~(np.isinf(np.abs(df.log_pval))) & ~np.isnan(df.log_foldchange) & ~np.isnan(df.log_pval)]

    pval_threshold = 1.5
    fold_threshold = 1.5

    df_negfold_pval = df[(df.log_pval > pval_threshold)  & (df.log_foldchange < -fold_threshold) ]
    df_posfold_pval = df[(df.log_pval > pval_threshold)  & (df.log_foldchange > fold_threshold) ]
    df_above_threshold = df[((df.log_pval <= pval_threshold)) | ((np.abs(df.log_foldchange) < fold_threshold)  )]

    out_dict = {}
    negfold_pval = {'substrate': list(df_negfold_pval.Substrate), 'foldchange': list(df_negfold_pval.log_foldchange), 'pval': list(df_negfold_pval.log_pval)}
    posfold_pval = {'substrate': list(df_posfold_pval.Substrate), 'foldchange': list(df_posfold_pval.log_foldchange), 'pval': list(df_posfold_pval.log_pval)}
    above_pval = {'substrate': list(df_above_threshold.Substrate), 'foldchange': list(df_above_threshold.log_foldchange), 'pval': list(df_above_threshold.log_pval)}
    out_dict['negfold_pval'] = negfold_pval
    out_dict['posfold_pval'] = posfold_pval
    out_dict['above_pval'] = above_pval
    os.remove(file_path)
    return out_dict

@app.route('/phosphoproteomics', methods = ['GET','POST'])
def phosphoproteomics():
    context = {}
    if request.method == 'POST':
        tmp_file_name = session['tmp_upload_file'] # get something from request
        session['tmp_upload_file'] = ''
        tmp_file_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_file_name)
        print(tmp_file_path)
        results_volcano = volcano(tmp_file_path)
        context['volcano'] = results_volcano
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



@app.route('/kinase/gene/<kin_name>.fasta')
def fasta_gene(kin_name):
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
