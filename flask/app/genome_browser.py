from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
from flask_wtf import CsrfProtect
import os
import pandas as pd
import queries




genome_browser_blueprint = Blueprint(
    'genomic_browser', __name__,
    template_folder = 'templates/genome_browser',
    url_prefix = '/phosphosites'
)


@genome_browser_blueprint.record
def record_params(setup_state):
  app = setup_state.app
  genome_browser_blueprint.config = dict([(key,value) for (key,value) in app.config.items()])


@genome_browser_blueprint.route('/<chromosome>', methods=['GET', 'POST'])
def genome_viewer_chrom(chromosome):
    #get requested values or assigned a default value
    requested_chrom_number =  int(chromosome)
    requested_chrom_ncbi = request.args.get('ncbi','NC_000001.11')
    requested_genom_start = int(request.args.get('gnm_start',104770758))
    requested_genom_end = int(request.args.get('gnm_end',104770760))

    phosphosites = queries.select_gral(genome_browser_blueprint.config['DATABASE'],'basic.reverse, upho.*', \
                        'uniprot_phosphosites upho LEFT JOIN basic_info basic ON basic.uniprot_id = upho.uniprot_id',\
                        'basic.chromosome LIKE "{0}" AND upho.genom_begin >= {1} AND upho.genom_end <= {2}'.format(requested_chrom_number,requested_genom_start, requested_genom_end))

    genom_browser_markers_list = []
    color = '0040FF' #blue

    for i in range(phosphosites.shape[0]):
        
        # check if the gene is in the reverse strand, and interchange start and end
        if (phosphosites.loc[i,'reverse']  == 'True'):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_end']) +':' + \
                                              str(phosphosites.loc[i,'genom_begin']) + '|' + \
                                              str(phosphosites.loc[i,'uniprot_id']) + '_'+\
                                              str(phosphosites.loc[i,'residue_position']) +' '+\
                                              '|' + color) 
        else:
            # IF not in the reverse (is in the forward), no need to interchange start and end
             genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_begin']) +':' + \
                                               str(phosphosites.loc[i,'genom_end']) + '|' + \
                                               str(phosphosites.loc[i,'uniprot_id']) + '_'+\
                                               str(phosphosites.loc[i,'residue_position']) + '|' + color)

    genom_browser_v = str(phosphosites['genom_begin'].min()-20) + ':' + str(phosphosites['genom_end'].max() +20)

    # string of markers
    genom_browser_markers = ','.join(genom_browser_markers_list)

        # generate url for ncbi genome viewer
    genom_browser_src ="?embedded=true&id="+requested_chrom_ncbi+"&tracks=[key:sequence_track,name:T16507,display_name:Sequence,id:T16507,dbname:SADB,annots:NA000001672.2,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:six_frames_translation,name:T11044,display_name:Six-frame translations,id:T11044,dbname:GenBank,annots:Six-frame translation,ShowOption:All,OrfThreshold:20,HighlightCodons:true,AltStart:false,shown:true,order:21][key:gene_model_track,name:T2000262,display_name:Genes\, NCBI Homo sapiens Annotation Release 105.20190906,id:T2000262,dbname:SADB,annots:NA000229419.1,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:22]&assm_context=GCF_000001405.13&mk="+genom_browser_markers+"&v="+genom_browser_v+"&c=ffff99&select=gi|224589800-0017f72a-001844c0-010a-ff7f5665-NA000229419.1;&slim=0&appname=pkinases"

    context = {'phosphosites':phosphosites, 'chromosome_number':requested_chrom_number,
               'requested_chrom_ncbi':requested_chrom_ncbi,'region_begin':requested_genom_start,
               'region_end':requested_genom_end, 'genom_browser_src':genom_browser_src}
    return render_template('/genome_browser.html',context=context)








@genome_browser_blueprint.route('/visualize/<chromosome>')
def genome_viewer_ncbi(chromosome):
    #get url arguments
    requested_chrom_number =  chromosome
    requested_chrom_ncbi = request.args.get('ncbi','NC_000001.11')
    requested_genom_start = int(request.args.get('gnm_start',104770758))
    requested_genom_end = int(request.args.get('gnm_end',104770760))

    phosphosites = queries.select_gral(genome_browser_blueprint.config['DATABASE'],'basic.reverse, upho.*', \
                        'uniprot_phosphosites upho LEFT JOIN basic_info basic ON basic.uniprot_id = upho.uniprot_id',\
                        'basic.chromosome LIKE "{0}" AND upho.genom_begin >= {1} AND upho.genom_end <={2}'.format(requested_chrom_number,requested_genom_start, requested_genom_end))


    genom_browser_markers_list = []
    color = '0040FF' #blue

    for i in range(phosphosites.shape[0]):
        
        # check if the gene is in the reverse strand, and interchange start and end
        if (phosphosites.loc[i,'reverse']  == 'True'):
            genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_end']) +':' + \
                                              str(phosphosites.loc[i,'genom_begin']) + '|' + \
                                              str(phosphosites.loc[i,'uniprot_id']) + '_'+\
                                              str(phosphosites.loc[i,'residue_position']) +' '+\
                                              '|' + color) 
        else:
            # IF not in the reverse (is in the forward), no need to interchange start and end
             genom_browser_markers_list.append(str(phosphosites.loc[i,'genom_begin']) +':' + \
                                               str(phosphosites.loc[i,'genom_end']) + '|' + \
                                               str(phosphosites.loc[i,'uniprot_id']) + '_'+\
                                               str(phosphosites.loc[i,'residue_position']) + '|' + color)

    genom_browser_v = str(phosphosites['genom_begin'].min()-20) + ':' + str(phosphosites['genom_end'].max() +20)

    # string of markers
    genom_browser_markers = ','.join(genom_browser_markers_list)

        # generate url for ncbi genome viewer
    genom_browser_src ="?embedded=true&id="+requested_chrom_ncbi+"&tracks=[key:sequence_track,name:T16507,display_name:Sequence,id:T16507,dbname:SADB,annots:NA000001672.2,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:six_frames_translation,name:T11044,display_name:Six-frame translations,id:T11044,dbname:GenBank,annots:Six-frame translation,ShowOption:All,OrfThreshold:20,HighlightCodons:true,AltStart:false,shown:true,order:21][key:gene_model_track,name:T2000262,display_name:Genes\, NCBI Homo sapiens Annotation Release 105.20190906,id:T2000262,dbname:SADB,annots:NA000229419.1,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:22]&assm_context=GCF_000001405.13&mk="+genom_browser_markers+"&v="+genom_browser_v+"&c=ffff99&select=gi|224589800-0017f72a-001844c0-010a-ff7f5665-NA000229419.1;&slim=0&appname=pkinases"
    

    context = {'genom_browser_src':genom_browser_src}
    return render_template('genome_chromosome_viewer.html', context=context)


