from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
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
import divide_sequences
import add_pubmed_link



phosphoproteomics_blueprint = Blueprint(
    'phosphoproteomics', __name__,
    template_folder = 'templates/phosphoproteomics',
    url_prefix = '/phosphoproteomics'
)

@phosphoproteomics_blueprint.record
def record_params(setup_state):
  app = setup_state.app
  phosphoproteomics_blueprint.config = dict([(key,value) for (key,value) in app.config.items()])


@phosphoproteomics_blueprint.route('/', methods = ['GET','POST'])
def phosphoproteomics():
    context = {}
    if request.method == 'POST':
        tmp_file_name = session['tmp_upload_file'] # get name of the file
        session['tmp_upload_file'] = ''
        fold_threshold = request.values['fc'] # get fold change threshold
        pval_threshold = request.values['pv'] # get pvalue change threshold
        tmp_file_path = os.path.join(phosphoproteomics_blueprint.config['UPLOAD_FOLDER'], tmp_file_name)

        try:
            output = phosphoproteomics_script.change_column_names(tmp_file_path)
            inhibitor = output['inhibitor']
            ddf = output['df']
            
            results_volcano = phosphoproteomics_script.volcano(ddf,pval_threshold, fold_threshold )
            df_volcano = phosphoproteomics_script.extract_above_threshold(ddf, results_volcano)
        
            query = "SELECT {} FROM {}".format('kinase, sub_gene, sub_mod_rsd, substrate', 'kinase_substrate')
            db = sqlite3.connect(phosphoproteomics_blueprint.config['DATABASE'])
            kin_substrate = pd.read_sql_query(query, db)
            db.close()
            z_score = phosphoproteomics_script.KSEA(ddf, kin_substrate) # for all
            z_score_volcano = phosphoproteomics_script.KSEA(df_volcano, kin_substrate) # just for the ones that are above threshold in volcano plot

            context['inhibitor'] = inhibitor.replace('_','').upper()
            context['volcano'] = results_volcano
            context['fold_threshold'] = fold_threshold
            context['pval_threshold'] = pval_threshold
            context['non_identified'] = z_score['non_identified']
            context['z_score'] = z_score['score']
            context['non_identified_volcano'] = z_score_volcano['non_identified']
            context['z_score_volcano'] = z_score_volcano['score']

        except:
            os.remove(tmp_file_path)
            return 'Impossible to calculate, something wrong in the input values. <a href="/"> Go back </a>'


    return render_template('phosphoproteomics.html', context = context)
