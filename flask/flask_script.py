# lines to change signaled with "modify:   !!"

from flask import Flask, render_template, redirect, request, url_for, session, flash
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


import pathlib
import re
path = str(pathlib.Path(__file__).parent.absolute())
DATABASE = re.sub(r'flask$','csv_tables/kinase_project.db',path)

UPLOAD_FOLDER = re.sub(r'flask$','csv_tables/',path)

app = Flask(__name__)
app.secret_key = 'my_secret_key'
csrf = CsrfProtect(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER



#### home
@app.route('/', methods=['GET','POST'])
def home():
    search_form = forms.Search_string()
    uploadfile_form = forms.UploadForm()
    phosphosite_form = forms.Phosphosite()
    inhibitor_form = forms.Inhibitors()
    context = {'uploadfile_form': uploadfile_form, 
               'search_form': search_form, 'phosphosite_form': phosphosite_form,
               'inhibitor_form':inhibitor_form}

    # kinase search form
    if request.method == 'POST' and search_form.search_string.data and search_form.validate_on_submit():
        requested_name = request.values['search_string']
        requested_name = requested_name.rstrip()
        if queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info', \
                                   'uniprot_id LIKE "{0}" OR name_human LIKE "{0}" OR prot_name LIKE "{0}"'.format(requested_name)):
            uniprot_id = queries.select_gral(DATABASE, 'uniprot_id','kinase_info',\
                                             'uniprot_id LIKE "{0}" OR name_human LIKE "{0}" OR prot_name LIKE "{0}"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase_data',kin_name= uniprot_id))
        elif queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info', \
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name)): # modify: get list of kinases
            uniprot_id = queries.select_gral(DATABASE, 'uniprot_id','kinase_info', \
                                             'uniprot_id LIKE "%{0}%"'.format(requested_name)).loc[0,'uniprot_id']
            return redirect(url_for('kinase_data',kin_name= uniprot_id))
        elif queries.query_n_results(DATABASE, 'uniprot_id','kinase_info',\
                                     'uniprot_id LIKE "%{0}%"'.format(requested_name))>1: # modify: get list of kinases
            return redirect(url_for('.kinase_search_result', search=requested_name))
        else:
            text_flash = "'"+ str(requested_name) + "'" + " not in our the database"
            flash(text_flash, 'kinase')
            return render_template('home.html', context = context)
        # modify:
        #  - if name is equal to a uniprot identifier redirect to /kinase/uniprotid
        # - else if requested_name correspond to one gene, redirect to the kinase/uniprotid
        # - else if requested_name one name nor to one gene, make a less restrictive querry and redirect to kinase_search_results


    # inhibitors form
    if request.method == 'POST' and inhibitor_form.data and inhibitor_form.validate_on_submit():
        inhib_requested = request.values['inhibitor_name']
        inhib_requested = inhib_requested.rstrip()
        if queries.query_is_unique(DATABASE, 'inn_name','inhibitors_gral_info', \
                                   'inn_name LIKE "{0}"'.format(inhib_requested)):
            inn_name = queries.select_gral(DATABASE, 'inn_name','inhibitors_gral_info',\
                                             'inn_name LIKE "{0}"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor_data',inhib_name= inn_name))
        elif queries.query_is_unique(DATABASE, 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "{0}"'.format(inhib_requested)): # modify: get list of kinases
            inn_name = queries.select_gral(DATABASE, 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "{0}"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor_data',inhib_name= inn_name))
        if queries.query_is_unique(DATABASE, 'inn_name','inhibitors_gral_info', \
                                   'inn_name LIKE "%{0}%"'.format(inhib_requested)):
            inn_name = queries.select_gral(DATABASE, 'inn_name','inhibitors_gral_info',\
                                             'inn_name LIKE "%{0}%"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor_data',inhib_name= inn_name))
        elif queries.query_is_unique(DATABASE, 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested)): # modify: get list of kinases
            inn_name = queries.select_gral(DATABASE, 'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested)).loc[0,'inn_name']
            return redirect(url_for('inhibitor_data',inhib_name= inn_name))
        elif queries.query_n_results(DATABASE,  'inn_name','inhibitors_gral_info',\
                                     'inn_name LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('.inhibitor_search_result', search=inhib_requested,type='inn'))
        elif queries.query_n_results(DATABASE,  'inn_name','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(inhib_requested))>1: # modify: get list of kinases
            return redirect(url_for('.inhibitor_search_result', search=inhib_requested,type='syn'))

        else:
            text_flash = "'"+ str(inhib_requested) + "'" + " not in our the database"
            flash(text_flash, 'inhibitor')
            return render_template('home.html', context = context)


    # uploaded file form
    if request.method == 'POST' and uploadfile_form.uploaded_file.data and uploadfile_form.validate_on_submit():
        file = request.files['uploaded_file']
        p_val_threshold = request.values['threshold_pval']
        threshold_foldchange = request.values['threshold_foldchange']
        random_name = str(random()).split('.')[1] #random number
        filename = random_name + '.tsv'
        session['tmp_upload_file'] = filename
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)

        return redirect(url_for('phosphoproteomics',fc=threshold_foldchange,pv=p_val_threshold), code = 307)

    
    return render_template('home.html', context = context)



#### kinases


@app.route('/kinase/results')
def kinase_search_result():
    # method: get , add filter
    requested_name = request.values['search']
    search_results = queries.select_gral(DATABASE, 'uniprot_id, name_human, chromosome, fasd_name, ensembl_gene_id','kinase_info', 'uniprot_id LIKE "%{0}%"'.format(requested_name))
    context = {'search_results':search_results}
    return render_template('kinase_search_results.html', context = context)


@app.route('/kinase/<kin_name>')
def kinase_data(kin_name):
    # modify: method: post, add filter !!
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info', "uniprot_id LIKE '{}'".format(kin_name)) : # modify: get list of kinases
        kin_name = kin_name
        gral_info = queries.select_gral(DATABASE, '*', 'kinase_info', 'uniprot_id LIKE "{}"'.format(kin_name))
        isoforms = list(queries.select_gral(DATABASE, 'isoform', 'isoforms', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'isoform'])

      # inhibitors
        inhibitors = list(queries.select_gral(DATABASE,'inn_name' ,'inhibitors_targets', 'targets LIKE "{}"'.format(gral_info.loc[0,'prot_name']) ).loc[:,'inn_name'])


        function_list_tmp = list(queries.select_gral(DATABASE, 'prot_function', 'kin_function', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'prot_function'])
        function_list =[]
        for function in function_list_tmp:
            function_list.append(add_pubmed_link.pubmed_link(function))

        # modify: add layers of information
        reactions_list = list(queries.select_gral(DATABASE, 'reaction_text', 'reactions', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
        cell_loc_list= list(queries.select_gral(DATABASE, 'DISTINCT subcell_location', 'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])

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
                   'targets':targets, 'phosphosites':phosphosites,
                   'inhibitors': inhibitors
                   }
        return render_template('kinase_data.html', context = context)
    else:
        return 'not found'


########################
# genome viewer
#

@app.route('/genome_viewer/<uniprot_id>')
def genome_viewer(uniprot_id):
    phosphosites = queries.select_gral(DATABASE, 'residue_position, modif, type_modif, genom_begin, genom_end', 'phosphosites', 'uniprot_id LIKE "{}"'.format(uniprot_id))

    gral_info = queries.select_gral(DATABASE, 'uniprot_id,reverse,chromosome', 'kinase_info', 'uniprot_id LIKE "{}"'.format(uniprot_id))
    chromosome_ncbi = queries.select_gral(DATABASE, 'ncbi_id','ncbi_chrom_id','chr LIKE  "{}"'.format(gral_info.loc[0,'chromosome'])).loc[0,'ncbi_id']

    genom_browser_markers_list = []
    color = '0040FF'
    genom_browser_v = ''
    
    if (gral_info.loc[0,'reverse']  == 'True'):
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_end']) +':' + \
            str(phosphosites.loc[i,'genom_begin']) + '|' + str(phosphosites.loc[i,'residue_position']) +' '+\
            phosphosites.loc[i,'type_modif']+ '|' + color) 
        genom_browser_v = str(phosphosites['genom_end'].min()-20)  + ':' +  str(phosphosites['genom_begin'].max() +20)
    else:
        color = '0040FF'
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_begin']) +':' + \
            str(phosphosites.loc[i,'genom_end']) + '|' + str(phosphosites.loc[i,'residue_position']) + '|' + color)
        genom_browser_v = str(phosphosites['genom_begin'].min()-20) + ':' + str(phosphosites['genom_end'].max() +20)

    genom_browser_markers = ','.join(genom_browser_markers_list)


    genom_browser_src ="?embedded=true&id="+chromosome_ncbi+"&tracks=[key:sequence_track,name:T16507,display_name:Sequence,id:T16507,dbname:SADB,annots:NA000001672.2,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:six_frames_translation,name:T11044,display_name:Six-frame translations,id:T11044,dbname:GenBank,annots:Six-frame translation,ShowOption:All,OrfThreshold:20,HighlightCodons:true,AltStart:false,shown:true,order:21][key:gene_model_track,name:T2000262,display_name:Genes\, NCBI Homo sapiens Annotation Release 105.20190906,id:T2000262,dbname:SADB,annots:NA000229419.1,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:22]&assm_context=GCF_000001405.13&mk="+genom_browser_markers+"&v="+genom_browser_v+"&c=ffff99&select=gi|224589800-0017f72a-001844c0-010a-ff7f5665-NA000229419.1;&slim=0&appname=pkinases"


    context = {'genom_browser_src':genom_browser_src}
    return render_template('genome_viewer.html', context=context)




@app.route('/genome_viewer/<chromosome>')
def genome_viewer_chrom(chromosome):
    phosphosites = queries.select_gral(DATABASE, 'residue_position, modif, type_modif, genom_begin, genom_end', 'phosphosites', 'uniprot_id LIKE "{}"'.format(uniprot_id))

    gral_info = queries.select_gral(DATABASE, 'uniprot_id,reverse,chromosome', 'kinase_info', 'uniprot_id LIKE "{}"'.format(uniprot_id))
    chromosome_ncbi = queries.select_gral(DATABASE, 'ncbi_id','ncbi_chrom_id','chr LIKE  "{}"'.format(gral_info.loc[0,'chromosome'])).loc[0,'ncbi_id']

    genom_browser_markers_list = []
    color = '0040FF'
    genom_browser_v = ''

    if (gral_info.loc[0,'reverse']  == 'True'):
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_end']) +':' + \
            str(phosphosites.loc[i,'genom_begin']) + '|' + str(phosphosites.loc[i,'residue_position']) +' '+\
            phosphosites.loc[i,'type_modif']+ '|' + color)
        genom_browser_v = str(phosphosites['genom_end'].min()-20)  + ':' +  str(phosphosites['genom_begin'].max() +20)
    else:
        color = '0040FF'
        for i in range(phosphosites.shape[0]):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_begin']) +':' + \
            str(phosphosites.loc[i,'genom_end']) + '|' + str(phosphosites.loc[i,'residue_position']) + '|' + color)
        genom_browser_v = str(phosphosites['genom_begin'].min()-20) + ':' + str(phosphosites['genom_end'].max() +20)

    genom_browser_markers = ','.join(genom_browser_markers_list)


    genom_browser_src ="?embedded=true&id="+chromosome_ncbi+"&tracks=[key:sequence_track,name:T16507,display_name:Sequence,id:T16507,dbname:SADB,annots:NA000001672.2,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:six_frames_translation,name:T11044,display_name:Six-frame translations,id:T11044,dbname:GenBank,annots:Six-frame translation,ShowOption:All,OrfThreshold:20,HighlightCodons:true,AltStart:false,shown:true,order:21][key:gene_model_track,name:T2000262,display_name:Genes\, NCBI Homo sapiens Annotation Release 105.20190906,id:T2000262,dbname:SADB,annots:NA000229419.1,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:22]&assm_context=GCF_000001405.13&mk="+genom_browser_markers+"&v="+genom_browser_v+"&c=ffff99&select=gi|224589800-0017f72a-001844c0-010a-ff7f5665-NA000229419.1;&slim=0&appname=pkinases"


    context = {'genom_browser_src':genom_browser_src}
    return render_template('genome_viewer.html', context=context)




#### inhibitors


@app.route('/inhibitor/search')
def inhibitor_search_result():
    # method: get, add filter
    requested_name = request.values['search']
    requested_type = request.values['type']
    if (requested_type == 'inn'):
        search_results = queries.select_gral(DATABASE,  'inn_name, phase, mw','inhibitors_gral_info',\
                                     'inn_name LIKE "%{0}%"'.format(requested_name))
    elif (requested_type == 'syn'):
        search_results = queries.select_gral(DATABASE,  'inn_name, phase, mw','inhibitors_synonims', \
                                     'synonyms LIKE "%{0}%"'.format(requested_name))
    context = {'search_results': search_results}
    return  render_template('inhibitor_search_results.html', context = context)

@app.route('/inhibitor/<inhib_name>')
def inhibitor_data(inhib_name):
    # method: post, add filter
    gral_info = queries.select_gral(DATABASE, '*','inhibitors_gral_info',\
                                    ' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = queries.select_gral(DATABASE, 'inh.targets, kin.uniprot_id',\
                                  'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                                  ' inn_name LIKE "{}"'.format(inhib_name))


    synonyms = list(queries.select_gral(DATABASE, 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    pdbid = list(queries.select_gral(DATABASE, 'pdbid','inhibitors_pdbid',\
                                      ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
    families = list(queries.select_gral(DATABASE, 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
    context = {'gral_info':gral_info, 'targets': targets, 'synonyms':synonyms, 'pdbid':pdbid, 'families':families}
    return render_template('inhibitor_data.html', context=context)



@app.route('/inhibitor/basic_info/<inhib_name>')
def inhibitor_data_information_iframe(inhib_name):
    # method: post, add filter
    gral_info = queries.select_gral(DATABASE, '*','inhibitors_gral_info',' inn_name LIKE "{}"'.format(inhib_name)).loc[0,:]
    targets = queries.select_gral(DATABASE, 'inh.targets, kin.uniprot_id','inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',' inn_name LIKE "{}"'.format(inhib_name))


    synonyms = list(queries.select_gral(DATABASE, 'synonyms','inhibitors_synonims',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
    pdbid = list(queries.select_gral(DATABASE, 'pdbid','inhibitors_pdbid',\
                                      ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
    families = list(queries.select_gral(DATABASE, 'kinase_families','inhibitors_kin_family',\
                                         ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
    context = {'gral_info':gral_info, 'targets': targets, 'synonyms':synonyms, 'pdbid':pdbid, 'families':families}
    return render_template('inhibitor_data_information.html', context=context)




### phosphoproteomics


@app.route('/phosphoproteomics', methods = ['GET','POST'])
def phosphoproteomics():
    context = {}
    if request.method == 'POST':
        tmp_file_name = session['tmp_upload_file'] # get name of the file
        session['tmp_upload_file'] = ''
        fold_threshold = request.values['fc'] # get fold change threshold
        pval_threshold = request.values['pv'] # get pvalue change threshold
        tmp_file_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_file_name)

        try:
            output = phosphoproteomics_script.change_column_names(tmp_file_path)
            inhibitor = output['inhibitor']
            ddf = output['df']
            
            results_volcano = phosphoproteomics_script.volcano(ddf,pval_threshold, fold_threshold )
            df_volcano = phosphoproteomics_script.extract_above_threshold(ddf, results_volcano)
        
            query = "SELECT {} FROM {}".format('kinase, sub_gene, sub_mod_rsd, substrate', 'kinase_substrate')
            db = sqlite3.connect(DATABASE)
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
    if queries.query_is_unique(DATABASE, 'uniprot_id', 'kinase_info','uniprot_id  LIKE "{}"'.format( kin_name)) : # modify: get list of kinases
        text = queries.select_gral(DATABASE, 'chromosome, genome_sequence, reverse, ensembl_gene_id, genome_starts, genome_ends','kinase_info', 'uniprot_id  LIKE "{}"'.format( kin_name))
        sequence = text.loc[0,'genome_sequence']
        divide_each = 40  # modify: change size !!
        seq_size = len(sequence)
        list_range = range(0,seq_size,divide_each)
        tmp_text= ''
        for i in list_range:
            tmp_text += sequence[i:i+divide_each] + '\n'
            seq_out = tmp_text.rstrip()
            header = '> '+ kin_name + ', length: ' + str(seq_size) + ', Chrom: ' + \
            text.loc[0,'chromosome'] + ', Reverse: ' + text.loc[0,'reverse'] + \
            ', Ensembl ID: ' + str( text.loc[0,'ensembl_gene_id']) + ', Start: ' + str(text.loc[0,'genome_starts']) + \
            ', Ends: ' + str(text.loc[0,'genome_ends'])
            text_out = '\n'.join([header, seq_out])
        return text_out, 200, {'Content-Type': 'text/plain; charset=utf-8'}
    else:
        return '', 200, {'Content-Type': 'text/plain; charset=utf-8'}


# this returns json !!
@app.route('/kinase/<kin_name>.json', methods=['GET'])
def kinase_json(kin_name):
    requested_columns = request.args.get('columns','')
    requested_columns = requested_columns.rstrip().lower()
    if queries.query_is_unique(DATABASE, 'uniprot_id','kinase_info',\
                                   'uniprot_id LIKE "{0}"'.format(kin_name)):
        requested_columns_list = requested_columns.split(',')
        if (('all' in requested_columns) or (requested_columns == '')):
            gral_info = queries.select_gral(DATABASE, '*','kinase_info', \
                                            'uniprot_id LIKE "{0}" '.format(kin_name))
            isoforms = queries.select_gral(DATABASE, 'kin.isoform, iso.prot_sequence',\
                                           'isoforms kin LEFT JOIN isoforms_info iso ON kin.isoform = iso.uniprot_id', \
                                           'kin.uniprot LIKE "{}"'.format(kin_name)) 
            inhibitors = queries.select_gral(DATABASE,\
                                                  'inh.inn_name, inhinfo.phase, inhinfo.mw, inhinfo.canonical_smiles, inhinfo.inchikey' ,\
                                                  'kinase_info kin  INNER JOIN inhibitors_targets inh ON kin.prot_name = inh.targets INNER JOIN inhibitors_gral_info inhinfo ON inh.inn_name = inhinfo.inn_name',\
                                                  'kin.uniprot_id LIKE "{}"'.format(kin_name))
            reactions_list = list(queries.select_gral(DATABASE, 'reaction_text', 'reactions',\
                                                      'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
            cell_loc_list= list(queries.select_gral(DATABASE, 'DISTINCT subcell_location',\
                                                    'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])
            diseases = queries.select_gral(DATABASE,\
                                           'DISTINCT disease_name, effect_text, disease_description', 'diseases',\
                                           'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
            targets = queries.select_gral(DATABASE, 'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa',\
                                          'kinase_substrate', 'kin_acc_id LIKE "{}"'.format(kin_name))
            phosphosites = queries.select_gral(DATABASE, 'residue_position, modif, type_modif, genom_begin, genom_end', \
                                               'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))
        
            output_dict = {'general_info':gral_info.to_dict('index')[0],
                           'isoforms': isoforms.to_dict('index'), 'inhibitors': inhibitors.to_dict('index'),
                           'reactions':reactions_list, 'cell_loc': cell_loc_list,
                           'diseases':diseases.to_dict('index'), 'targets':targets.to_dict('index'),
                           'phosphosites': phosphosites.to_dict('index')}
        else:
            output_dict = {}
            for column  in requested_columns_list:
                if ('info' in column):
                    gral_info = queries.select_gral(DATABASE, '*','kinase_info', \
                                                    'uniprot_id LIKE "{0}" '.format(kin_name))
                    output_dict['general_info'] = gral_info.to_dict('index')[0]
                elif ('isoforms' in column):
                    isoforms = queries.select_gral(DATABASE, 'kin.isoform, iso.prot_sequence',\
                                                   'isoforms kin LEFT JOIN isoforms_info iso ON kin.isoform = iso.uniprot_id', \
                                                   'kin.uniprot LIKE "{}"'.format(kin_name)) 
                    output_dict['isoforms'] = isoforms.to_dict('index')
                elif ('inhibitors' in column):
                    inhibitors = queries.select_gral(DATABASE,\
                                                     'inh.inn_name, inhinfo.phase, inhinfo.mw, inhinfo.canonical_smiles, inhinfo.inchikey' ,\
                                                     'kinase_info kin  INNER JOIN inhibitors_targets inh ON kin.prot_name = inh.targets INNER JOIN inhibitors_gral_info inhinfo ON inh.inn_name = inhinfo.inn_name',\
                                                     'kin.uniprot_id LIKE "{}"'.format(kin_name))
                    output_dict['inhibitors'] = inhibitors.to_dict('index')
                elif ('reactions' in column):
                    reactions_list = list(queries.select_gral(DATABASE, 'reaction_text', 'reactions',\
                                                              'uniprot LIKE "{}"'.format(kin_name)).loc[:,'reaction_text'])
                    output_dict['reactions'] = reactions_list
                elif ('cell_loc' in column):
                    cell_loc_list = list(queries.select_gral(DATABASE, 'DISTINCT subcell_location',\
                                                             'subcell_location', 'uniprot LIKE "{}"'.format(kin_name)).loc[:,'subcell_location'])
                    output_dict['cell_loc'] = cell_loc_list
                elif ('diseases' in column):
                    diseases = queries.select_gral(DATABASE,\
                                                   'DISTINCT disease_name, effect_text, disease_description', 'diseases',\
                                                   'uniprot LIKE "{}" AND disease_name NOT LIKE "" ORDER BY disease_name'.format(kin_name))
                    output_dict['diseases'] = diseases.to_dict('index')
                elif ('phosphosites_self' in column):
                    phosphosites = queries.select_gral(DATABASE, 'residue_position, modif, type_modif, genom_begin, genom_end', \
                                                       'phosphosites', 'uniprot_id LIKE "{}"'.format(kin_name))
                    output_dict['phosphosites_self'] = phosphosites.to_dict('index')
                elif ('phosphosites_targets' in column):
                    targets = queries.select_gral(DATABASE, 'sub_acc_id, sub_gene, sub_mod_rsd, site_7_aa',\
                                                  'kinase_substrate', 'kin_acc_id LIKE "{}"'.format(kin_name))
                    output_dict['phosphosites_targets'] = targets.to_dict('index')
        return jsonify(output_dict)
    

    else:
        return ''

## inhibitors json api
@app.route('/inhibitor/<inhib_name>.json', methods=['GET'])
def inhibitor_json(inhib_name):
    requested_columns = request.args.get('columns','')
    requested_columns = requested_columns.rstrip().lower()
    if queries.query_is_unique(DATABASE, 'inn_name','inhibitors_gral_info', \
                               'inn_name LIKE "{0}"'.format(inhib_name)):
        requested_columns_list = requested_columns.split(',')
        if (('all' in requested_columns) or (requested_columns == '')):

            gral_info = queries.select_gral(DATABASE, 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                            'inhibitors_gral_info',\
                                            ' inn_name LIKE "{}"'.format(inhib_name))
            targets = queries.select_gral(DATABASE, 'inh.targets, kin.uniprot_id',\
                                      'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                                      ' inn_name LIKE "{}"'.format(inhib_name))
            synonyms_list = list(queries.select_gral(DATABASE, 'synonyms','inhibitors_synonims',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
            pdbid_list = list(queries.select_gral(DATABASE, 'pdbid','inhibitors_pdbid',\
                                              ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
            families_list = list(queries.select_gral(DATABASE, 'kinase_families','inhibitors_kin_family',\
                                                 ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
        
            output_dict = {'general_info':gral_info.to_dict('index')[0],
                           'targets':targets.to_dict('index'),
                           'synonyms':synonyms_list, 'Protein_databank_id':pdbid_list,
                           'families':families_list
                            }
            return jsonify(output_dict)
        else:
            output_dict = {}
            for column  in requested_columns_list:
                if ('info' in column):
                    gral_info = queries.select_gral(DATABASE, 'inn_name, phase,mw, image_url, canonical_smiles, inchikey',\
                                                    'inhibitors_gral_info',\
                                                    ' inn_name LIKE "{}"'.format(inhib_name))
                    output_dict['general_info'] = gral_info.to_dict('index')[0]
                elif ('targets' in column):
                    targets = queries.select_gral(DATABASE, 'inh.targets, kin.uniprot_id',\
                          'inhibitors_targets inh LEFT JOIN kinase_info kin ON kin.prot_name = inh.targets',\
                          ' inn_name LIKE "{}"'.format(inhib_name))
                    output_dict['targets'] = targets.to_dict('index')
                elif ('synonyms' in column):
                    synonyms_list = list(queries.select_gral(DATABASE, 'synonyms','inhibitors_synonims',\
                                     ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'synonyms'])
                    output_dict['synonyms'] = synonyms_list
                elif ('pbid' in column):
                    pdbid_list = list(queries.select_gral(DATABASE, 'pdbid','inhibitors_pdbid',\
                                  ' inn_name LIKE "{}"'.format(inhib_name)).loc[:, 'pdbid'])
                    output_dict['Protein_databank_id'] = pdbid_list
                elif ('families' in column):
                    families_list = list(queries.select_gral(DATABASE, 'kinase_families','inhibitors_kin_family',\
                                     ' inn_name LIKE "{}"'.format(inhib_name)).loc[:,'kinase_families'])
                    output_dict['families'] = families_list
            return jsonify(output_dict)
    else:
        return ''

    

if __name__ == '__main__':
    app.run(debug=True, port = 8000)
