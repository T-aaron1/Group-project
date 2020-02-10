#!/usr/bin/env python
# coding: utf-8

# In[16]:


import pandas as pd
import numpy as np
import re
import os

#r'C:\Users\Timot\Downloads\mux.tsv'

file=('~/Escritorio/mux.tsv')

def file_handle(file):
    ''' reads phosphoproteomics data, returns a dictionary with two entries: substrate-position and
        another entry that contains the foldchange information for the inhibitors. Each inhibitor dataframe is
        inside a dictionary with the name of the inhibitor. This output is to be used as input of other functions
        for phosphoproteomics analysis.
    '''
    df = pd.DataFrame(pd.read_csv(file, sep='\t'))
    #os.remove(file_path)
    df.columns = map(str.lower, df.columns)
    #drop rows and columns with no value
    df.dropna(axis=1, how='all', inplace=True)
    df.dropna(axis=0, how='all', inplace=True)
    #split substrate column into name and position 
    df['subst_position'] = df.substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    df['subst_name'] = df.substrate.str.split("(", n=1, expand=True)[0]
    substrate = df[['subst_name','subst_position']]
    #extracting unique inhibitor name 
    inhibitors_list = []
    inhibitors_dict={}
    for name in df.columns:
        if ('fold_change' in name):
                inhibitor = name.rsplit('_')[0]
                inhibitor = inhibitor.lower().rstrip()
                inhibitors_list.append(inhibitor)
    for i in range(len(inhibitors_list)):
        tmp_inhib_df=df.filter(regex= inhibitors_list[i])
        inhibitor = inhibitors_list[i] + '_'
        tmp_inhib_df.columns = [inh.replace(inhibitor,'') for inh in tmp_inhib_df.columns]
        inhibitors_dict[inhibitors_list[i]]= tmp_inhib_df
    output ={'substrate': substrate, 'inhibitors_dict' :inhibitors_dict }
    return output
