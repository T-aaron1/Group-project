from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
#from flask_csv import send_csv
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


documentation_blueprint = Blueprint(
    'documentation', __name__,
    template_folder = 'templates/documentation',
    url_prefix = '/documentation'
)



### Documentation

@documentation_blueprint.route('/general')
def documentation_general():
    return render_template('documentation_general.html')

@documentation_blueprint.route('/api')
def documentation_api():
    return render_template('documentation_api.html')

@documentation_blueprint.route('/stats')
def documentation_stats():
    x = [1,2,3]
    return render_template('documentation_stats.html', x = x, y = x)

