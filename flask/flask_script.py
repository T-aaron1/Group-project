from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
from flask import Response # for api: fasta , csv and so on
from flask import jsonify
from flask_wtf import CsrfProtect
import os
import pandas as pd
import numpy as np
import sqlite3
from random import random
import queries
import divide_sequences
import add_pubmed_link
from app.home import main_blueprint
from app.kinases import kinase_blueprint
from app.inhibitors import inhibitor_blueprint
from app.phosphoproteomics import phosphoproteomics_blueprint
from app.documentation import documentation_blueprint
from app.genome_browser import genome_browser_blueprint


import pathlib
import re
path = str(pathlib.Path(__file__).parent.absolute())
DATABASE = re.sub(r'flask$','csv_tables/kinase_project.db',path)


# folder where uploaded files are going to be temporarily saved
UPLOAD_FOLDER = re.sub(r'flask$','csv_tables/',path)


app = Flask(__name__, static_folder = './static')
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['DATABASE'] = DATABASE


# registering bluprints
app.register_blueprint(main_blueprint)
app.register_blueprint(genome_browser_blueprint)
app.register_blueprint(kinase_blueprint)
app.register_blueprint(inhibitor_blueprint)
app.register_blueprint(phosphoproteomics_blueprint)
app.register_blueprint(documentation_blueprint)

if __name__ == '__main__':
    app.run(debug=True, port = 8000)
