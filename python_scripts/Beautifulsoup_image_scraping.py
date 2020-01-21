#!/usr/bin/env python
# coding: utf-8

# In[2]:


import requests
r = requests.get('http://www.icoa.fr/pkidb/')


# In[2]:


print (r.text)


# In[6]:


#parse using BeautifoulSoup
from bs4 import BeautifulSoup
soup = BeautifulSoup(r.text, 'html.parser')


# In[7]:


results = soup.find_all('tr')


# In[8]:


print(results)


# In[9]:


len(results)


# In[10]:


print(results[0])


# In[11]:



final_results = results[1:]
print(final_results[0])


# In[12]:


first_result =results[1]
print(first_result)


# In[13]:


first_result.find('td').text


# In[14]:


first_result.find('td').contents


# In[15]:


first_result.find('td').contents[2]['src']


# In[16]:


chemical_str = []
for i in final_results:
    name = i.find('td').text
    image = i.find('td').contents[2]['src']
    image = "www.icoa.fr/pkidb/" + image
    print(image)
    
    chemical_str.append((name, image))


# In[17]:


print(chemical_str)


# In[18]:


#applying a tapular data structure 
#let's convert the tupules into a df
import pandas as pd
df = pd.DataFrame(chemical_str, columns=['name', 'image'])


# In[19]:


df.head()


# In[20]:


new_df= df.sort_values(by='name')
print(new_df)


# In[21]:


new_df.to_csv("Chemical_str_inhibitor.csv")

