from flask import Flask, Blueprint, render_template


documentation_blueprint = Blueprint(
    'documentation', __name__,
    template_folder = 'templates/documentation',
    url_prefix = '/documentation'
)



### Documentation

@documentation_blueprint.route('/phosphoproteomics')
def documentation_general():
    return render_template('documentation_phosphoproteomics.html')

@documentation_blueprint.route('/api')
def documentation_api():
    return render_template('documentation_api.html')

