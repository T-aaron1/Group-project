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
        cv_treatment_threshold = request.values['cvt']
        tmp_file_path = os.path.join(phosphoproteomics_blueprint.config['UPLOAD_FOLDER'], tmp_file_name)

        try:
            # get the data from the database
            query = "SELECT {} FROM {}".format('basick.prot_name AS kinase, basicsub.gene AS sub_gene, subs.sub_mod_rsd, basicsub.prot_name AS substrate', \
                                               'kinase_substrate subs LEFT JOIN basic_info basick ON subs.kin_acc_id = basick.uniprot_id LEFT JOIN  basic_info basicsub ON subs.sub_acc_id = basicsub.uniprot_id')
            db = sqlite3.connect(phosphoproteomics_blueprint.config['DATABASE'])
            kin_substrate = pd.read_sql_query(query, db)
            db.close()

            output = phosphoproteomics_script.file_handle(tmp_file_path, kin_substrate)
            substrate  = output['substrate']

        except:
            os.remove(tmp_file_path)
            return 'Something went wrong during initial file handling, please check that the file has the right format and column names <a href="/"> Go back </a>'

        try :
            print(output['inhibitors_dict'].keys())
            print('------------------ññ----')
            results_dict = {}
            for key in output['inhibitors_dict'].keys():
                inhibitor = key
                ddf = substrate.join(output['inhibitors_dict'].get(key))
                results_volcano = phosphoproteomics_script.volcano(ddf,pval_threshold, fold_threshold,  cv_treatment_threshold)
                df_volcano = phosphoproteomics_script.extract_above_threshold(ddf, results_volcano)

                # just in case the it was empty
                if (('fold_change' in ddf.columns) and (ddf.shape[0]>0)):
                    z_score = phosphoproteomics_script.KSEA(ddf) # for all
                else:
                    z_score =  {'non_identified':[],'score':[]}
                    # just for the ones that are above threshold in volcano plot
                    # cause of the thresholds, sometimes there could be empty dataframes
                if (('fold_change' in df_volcano.columns) and (df_volcano.shape[0] > 0)):
                    z_score_volcano = phosphoproteomics_script.KSEA(df_volcano)
                else:
                    z_score_volcano = {'non_identified':[],'score':[]}

                non_identified = z_score['non_identified'] if z_score['non_identified'] else 0
                non_identified_volcano =  z_score_volcano['non_identified'] if z_score_volcano['non_identified'] else 0
                results_dict[inhibitor] = {
                    'volcano': results_volcano, 'fold_threshold': fold_threshold,
                    'pval_threshold': pval_threshold,
                    'non_identified': non_identified,
                    'z_score': z_score['score'],
                    'non_identified_volcano': non_identified_volcano,
                    'z_score_volcano': z_score_volcano['score']
                    }
                context['results'] = results_dict
                context['n_results'] = len(results_dict.keys())

            return render_template('phosphoproteomics.html', context = context)
        except:
            return 'Something went wrong during data analysis, please check that the file has the right format and column names <a href="/"> Go back </a>'
