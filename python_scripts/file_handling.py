#!/usr/bin/env python
# coding: utf-8

# In[16]:


import pandas as pd
import numpy as np
import re
import os

file=(r'C:\Users\Timot\Downloads\mux.tsv')

def file_handle(tmp_file):
    df = pd.DataFrame(pd.read_csv(tmp_file, sep='\t'))
    
    #os.remove(file_path)
    
    df.columns = map(str.lower, df.columns)
    
    #drop rows and columns with no value
    df.dropna(axis=1, how='all', inplace=True)
    df.dropna(axis=0, how='all', inplace=True)
    #substrate  control_mean          mean  fold_change   p-value    ctrlCV   treatCV



    #split substrate column into name and position 
    df['subst_position'] = df.substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    df['subst_name'] = df.substrate.str.split("(", n=1, expand=True)[0]

    #creating a new dataframe with control_mean, substrate position and substrate name


    #control_subst = df(['subst_name','subst_position'])
    #df.drop(['subst_name','subst_position','control_mean'], axis=1, inplace=True)



    #extracting unique inhibitor name 
    inhibitors = []
    for name in df.columns:
        if ('fold_change' in name):
                inhibitor = name.rsplit('_')[0]
                inhibitor = inhibitor.lower().rstrip()
                inhibitors.append(inhibitor)
    #print (inhibitors)
     
    inhibitors_dict={}
    for i in range(len(inhibitors)):
        X=df.filter(regex= inhibitors[i])
        
        inhibitors_dict[inhibitors[i]]= X
    
    df['subst_position'] = df.substrate.str.split("(", n=1, expand=True)[1].str.replace(")", "")
    df['subst_name'] = df.substrate.str.split("(", n=1, expand=True)[0]

    control_subst = df[['subst_name','subst_position','control_mean']]
 
    #needs to iteratie over the dictionary and change the column header value of each inhibitor in inhibitor_dict 
    for j in range(len(inhibitors_dict)):
        for name in df.columns:
            df.rename(columns = {name: name.lower()}, inplace = True)
        for name in df.columns:
            if inhibitors in name.lower():
                df.rename(columns = {name: name.lower().replace(inhibitors,'')}, inplace = True)
                df.dropna(axis= 1, how='all', inplace = True)
    
    
    data_dict ={'control_subst': control_subst, 'inhibitors_dict' :inhibitors_dict }
    
    
    
    
    
    
    
    return data_dict


# In[17]:


file_handle(file)


# In[ ]:





# In[ ]:




