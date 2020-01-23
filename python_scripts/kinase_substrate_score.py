#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


#substrates_file = "/homes/dtg30/Desktop/Kinase_Substrate_Dataset.csv"
substrates_file = "/homes/ta317/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/Kinase_Substrate_Dataset.csv"

data = pd.read_csv(substrates_file)
data.dropna(axis= 1, how='all', inplace = True)


#upload_file = "/homes/ta317/Desktop/az20.tsv"
upload_file = "/homes/ta317/Downloads/az20.tsv"

upload = pd.read_csv(upload_file, sep = '\t')


upload.dropna(axis= 1, how='all', inplace = True)


upload['subst_position'] = ''
upload['subst_position'] = upload.Substrate.str.split("(", n=1, expand = True)[1].str.replace(")","")
upload['subst_name'] = upload.Substrate.str.split("(", n=1, expand = True)[0]


#get kinase name for substrate / position

# possibly solve this with a join, instead of a for loop, mainly beacause of the time it takes


upload2 = upload
#upload2 = upload2.join(data[['KINASE','SUBSTRATE','SUB_MOD_RSD']].set_index(['SUBSTRATE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )

upload2 = upload2.join(data[['KINASE','SUB_GENE','SUB_MOD_RSD']].set_index(['SUB_GENE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )

upload2 = upload2.join(data[['KINASE','SUB_ACC_ID','SUB_MOD_RSD']].set_index(['SUB_ACC_ID','SUB_MOD_RSD']), on= ['subst_name','subst_position'] , rsuffix='_sub_acc')


upload2.columns
upload2['kinase_final'] = upload2[['KINASE','KINASE_sub_acc']].fillna('').sum(axis=1)


upload3=upload2.drop(upload2.ix[:, 'subst_name':'kinase_final'].columns, axis = 1) 
    
    #upload2.drop(['KINASE','KINASE_sub_acc'])


# In[2]:


#print (upload2)
#upload2.kinase_final[1:20]
upload3=upload2.drop(upload2.ix[:, 'KINASE':'KINASE_sub_acc'].columns, axis = 1) 


# In[3]:


print (upload3)


# In[4]:


upload3.to_csv(r'/homes/ta317/Desktop/group_proj/venv/src/Group-project/csv_tables/kinases/az20_kinase_substrate.csv')


# In[5]:


upload3.kinase_final.unique()


# In[6]:


import numpy as np

# Drop rows with null values
upload3 = upload3.dropna(axis=0)

#Drop rows with inf value
upload3 = upload3[upload3.AZ20_fold_change != 'inf']

#log2 the substrate to calculate score for kinase activity
upload3['Log2Substrate_fold_change'] = np.log2(upload3.AZ20_fold_change) 

#upload4 = upload3[upload3.AZ20_fold_change != "inf"]

print (upload3)
upload3.to_csv("test.csv", index=False)


# In[7]:


upload3['mean_log2(FC)']=upload3["Log2Substrate_fold_change"].mean()

print (upload3)


# In[9]:


#m denotes the total number of phosphosite substrates identified from the experiment 
#that annotate to the specified kinase

#upload3.count()
kinase_count =upload3['kinase_final'].value_counts()

#upload3['kinase_final'].replace('', np.nan, inplace=True)
#upload3.dropna(subset=['kinase_final'], inplace=True)

kinase_count =upload3['kinase_final'].value_counts()

#print (upload3)
kinase_count_df = pd.DataFrame(kinase_count) 


kinase_count_df['sqrt_kinase'] = np.sqrt(kinase_count_df) 

kinase_count_df

#upload3 = upload2.join(data[['KINASE','SUB_GENE','SUB_MOD_RSD']].set_index(['SUB_GENE','SUB_MOD_RSD']), on= ['subst_name','subst_position'] )
#kinase_count_df.to_csv("kinase_count_df.csv", index=False)


# In[10]:


upload3


# In[17]:


upload3.groupby('kinase_final')['Log2Substrate_fold_change'].mean()


# In[ ]:





# In[48]:


#kinase activity calculation table
kinase_count_df['Standard_deviation']= np.std(upload3.Log2Substrate_fold_change)


kinase_count_df['Log2Substrate_fold_change'] =upload3.groupby('kinase_final')['Log2Substrate_fold_change'].mean()

kinase_count_df['mean_log2(FC)'] = upload3['mean_log2(FC)']=upload3["Log2Substrate_fold_change"].mean()

kinase_count_df.drop(kinase_count_df.index[0])

kinase_count_df


# In[35]:


#Calculating the relative kinase activity score

kinase_count_df['KSEA'] =  (kinase_count_df['meanlogfc'] - kinase_count_df['mean_log2(FC)']*
                            kinase_count_df['sqrt_kinase'])/ kinase_count_df['Standard_deviation']

kinase_count_df

#kinase_count_df.to_csv("score.csv")


# In[52]:


from scipy.stats import norm

kinase_count_df['P_value']= norm.sf(abs(kinase_count_df['KSEA']))*2

kinase_count_df.sort_values(['P_value'], ascending = True)

test_plot=kinase_count_df[kinase_count_df['P_value']<0.05]

test_plot.to_csv("test_plot.csv")

kinase_count_df.to_csv("testing.csv")

