import flask
from flask_wtf import CsrfProtect
from flask import request, render_template,jsonify, redirect


app = flask.Flask(__name__)
app.config["DEBUG"] = True
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)


import forms

## form
@app.route('/form', methods = ['GET', 'POST'])
def formulario():
    form_username = forms.TmpForm(request.form)

    if request.method == 'POST' and form_username.validate():
        print(form_username.uploaded_file.data)
        return redirect('/resultado')
    context = {}
    context['form_username'] = form_username
    return render_template('tmp_home.html', context = context)


@app.route('/resultado',  methods = ['GET', 'POST'])
def resultado():
        return render_template('phosphoproteomics.html')

##########






# Create some test data for our catalog in the form of a list of dictionaries.
booksl = [
    {'id': 0,
     'title': 'A Fire Upon the Deep',
     'author': 'Vernor Vinge',
     'first_sentence': 'The coldsleep itself was dreamless.',
     'year_published': '1992'},
    {'id': 1,
     'title': 'The Ones Who Walk Away From Omelas',
     'author': 'Ursula K. Le Guin',
     'first_sentence': 'With a clamor of bells that set the swallows soaring, the Festival of Summer came to the city Omelas, bright-towered by the sea.',
     'published': '1973'},
    {'id': 2,
     'title': 'Dhalgren',
     'author': 'Samuel R. Delany',
     'first_sentence': 'to wound the autumnal city.',
     'published': '1975'}
]


@app.route('/<bla>.fasta', methods=['GET'])
def api(bla):
    text1 = '<p>> blkasdoiasj</p>'
    text2 = '<p>isadoiasj</p>'
    text_out = ''.join([text1,text2])
    return text_out


# !! this does the api thing with .fasta return
from flask import Response
@app.route('/ajax_ddl.fasta')
def ajax_ddl():
    x = "some data you\nwasasdasdoiasjhoijnt to"
    return x, 200, {'Content-Type': 'text/plain; charset=utf-8'} 








@app.route('/', methods=['GET'])
def home():
    return '''<h1>Distant Reading Archive</h1>
<p>A prototype API for distant reading of science fiction novels.</p>'''


# this returns json !!
@app.route('/api/v1/resources/books/all.json', methods=['GET'])
def api_all():
    books =     [{'id': 2,
     'title': 'Dhalgren',
     'author': 'Samuel R. Delany',
     'first_sentence': 'to wound the autumnal city.',
     'published': '1975'}]
    return jsonify(books)


from flask import Response
@app.route('/test.xml')
def test_xml():
#    xml = 'foo'
    r = Response(response="TEST OK", status=200, mimetype="application/xml")
    r.headers["Content-Type"] = "text/xml; charset=utf-8"
    return r
#    return Response(xml, mimetype='text/xml')




@app.route('/api/v1/resources/books.json', methods=['GET'])
def api_id():
    # Check if an ID was provided as part of the URL.
    # If ID is provided, assign it to a variable.
    # If no ID is provided, display an error in the browser.
    if 'id' in request.args:
        id = int(request.args['id'])
    else:
        return "Error: No id field provided. Please specify an id."

    # Create an empty list for our results
    results = []

    # Loop through the data and match results that fit the requested ID.
    # IDs are unique, but other fields might return many results
    for book in books:
        if book['id'] == id:
            results.append(book)

    # Use the jsonify function from Flask to convert our list of
    # Python dictionaries to the JSON format.
    return jsonify(results)

if __name__ == '__main__':
    app.run(debug = True, port = 8001)
