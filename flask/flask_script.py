# lines to change signaled with "modify:   !!"

from flask import Flask, render_template, redirect, request
import forms
#from flask_csv import send_csv
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect



app = Flask(__name__)
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)



#### home
@app.route('/', methods=['GET','POST'])
def home():
    uploadfile_form = forms.UploadPhosphoproteomics(request.form)
    context = {'uploadfile_form': uploadfile_form}

    if request.method == 'POST' and uploadfile_form.validate():
        print(uploadfile_form.uploaded_file.data)
        return redirect('/phosphoproteomics')
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
@app.route('/phosphoproteomics', methods = ['GET','POST'])
def phosphoproteomics():
    data_form = forms.UploadPhosphoproteomics(request.form)
    return render_template('phosphoproteomics.html')

### Documentation

@app.route('/documentation/general')
def documentation_general():
    return render_template('documentation_general.html')

@app.route('/documentation/api')
def documentation_api():
    return render_template('documentation_api.html')

@app.route('/documentation/stats')
def documentation_stats():
    return render_template('documentation_stats.html')





### API

# return fasta file
# !! this does the api thing with .fasta return
# this could be done also to generate a csv

@app.route('/kinase/<kin_name>.fasta')
def ajax_ddl(kin_name):
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
def ajax_ddl(kin_name):
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
