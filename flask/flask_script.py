from flask import Flask, render_template, redirect, request
import forms

app = Flask(__name__)

#### home
@app.route('/', methods=['GET','POST'])
def home():
    form =  forms.UploadPhosphoproteomics()
    context = {'form': form}

    data_form = forms.UploadPhosphoproteomics(request.form)

    if request.method == 'POST' and data_form.validate():
        print(data_form.uploaded_file.data)
        render_template('phosphoproteomics.html')
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
    # method: post, add filter
    kin_name = kin_name
    context = {'kin_name':kin_name}
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

    if request.method == 'POST' and data_form.validate():
        print(data_form.uploaded_file.data )
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




if __name__ == '__main__':
    app.run(debug=True, port = 8000)

