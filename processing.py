#!/usr/bin/env python
# coding: utf-8

# # Post processing scraped data
# 
# This notebook processes the scraped data from Portal de la Reserca to create the Nodelist and Edgelist to plot in Gephi.

# # Import modules

# In[16]:


import numpy as np
import pandas as pd
from ast import literal_eval


# # Get researchers of interest

# In[17]:


res_df = pd.read_csv('./data/nodelist.csv')
res_df_IGTP = res_df.loc[res_df['institution'] == 'IGTP']
res_IGTP = res_df_IGTP['id'].unique()


# # Get edgelist

# In[18]:


# papers0_df = pd.read_csv('./data/papers_0.csv')
# papers1_df = pd.read_csv('./data/papers_1.csv')
# papers_df = papers0_df.append(papers1_df)
# papers_df = papers_df.drop_duplicates()
# papers_df.to_csv('./data/papers.csv', index=False)


# In[19]:


papers_df = pd.read_csv('./data/papers.csv')


# In[22]:


# papers_df_backup = papers_df.copy()
papers_df = papers_df_backup.copy()


# # Merge

# In[23]:


# Convert strings to lists
def convert_to_list(x):
    try:
        result = literal_eval(x)
    except ValueError:
        result = np.nan
    return result

papers_df['orcids'] = papers_df['orcids'].apply(lambda x: convert_to_list(x))

# Check if any coauthor is in institution
def belongs_to_list(row, df):
    try:
        result = bool(set(row['orcids']) & set(df))
    except TypeError:
        result = False
    return result

mask = papers_df.apply(lambda x: belongs_to_list(x, res_IGTP), axis=1)

selected_df = papers_df[mask]


# # Todo
# 1. Create boolean vector P for each paper of len(authors) that identifies coauthors.
# 2. Create boolean vector C of coauthor focus. inner(C,P) = 1 if there is collaboration.

# ### Create bollean vector P for each paper

# In[135]:


# Test with n authors for debugging
# n = 5


# In[157]:


# Get papers of authors in instutition
papers = selected_df['orcids'].copy()


# In[158]:


# Create subselection for debugging
# papers = papers[:n]


# In[159]:


papers = papers.reset_index(drop=True)


# In[160]:


# Get unique list of authors
authors_index = list(set(papers.sum()))
authors_index.sort()


# In[161]:


# Create boolean matrix with papers
paper_bool_df = pd.DataFrame(index=authors_index)


# In[162]:


for i, paper in enumerate(papers):
    paper_bool_df.loc[:,i] = 0
    for orcid in paper:
        paper_bool_df.loc[orcid, i] = 1


# ### Create boolean matrix of coauthor combinations

# In[163]:


# Identify coauthors of the first paper
# paper = paper_bool_df[[0]].values.flatten()
# np.where(paper == 1)


# In[164]:


# Create custom vector to check
# C = np.zeros(len(authors_index))
# C[61] = 1
# C[97] = 1


# In[165]:


# Check inner product
# paper.dot(C)


# In[166]:


# M = np.zeros((3,9))


# In[167]:


# M[:,:3] = np.identity(3)


# In[ ]:


# n = 4
n = len(authors_index)

combinations_mat = np.zeros((n,n*n))

# Creation combination matrix
C0 = np.identity(n)
for i in range(n):
    print(f"Progress: {i/n*100:.0f}%.", end="\r")
    C = C0
    # Set rows
    C[i] = 1
    if i>=1:
        C[i-1] = 0
    # Set columns
        for j in range(i):
            C[:,j] = 0
        
    combinations_mat[:,i*n:(i+1)*n] = C
        
    C0 = C


# In[ ]:


combinations_mat


# In[ ]:


papers_mat = paper_bool_df.to_numpy()


# In[ ]:


result = np.dot(papers_mat, combinations_mat)


# In[ ]:


papers_mat.shape


# In[ ]:


combinations_mat.shape


# In[ ]:


result = result-1
result = result.clip(0)
links = result.sum(axis=1)


# In[ ]:


links


# In[1]:


# links


# In[2]:


# paper_bool_df


# In[3]:


# papers_test = np.array([[1,1,0,0],[1,1,0,0], [1,0,1,0]])


# In[4]:


# papers_test


# In[5]:


# result = np.dot(papers_test, M)
# result


# In[6]:


# result = result-1
# result


# In[7]:


# result = result.clip(0)
# result


# In[8]:


# result.sum(axis=0)


# In[9]:


# max(result-1,0)


# In[416]:


# C0.shape


# In[10]:


# C1 = np.identity(4)
# C[0] = 1

# C2 = np.identity(4)
# C2[1]=1
# C2[:,0]=0
# C2

# C3 = C2
# C3[1] = 0
# C3[2] = 1
# C3[:,0]=0
# C3[:,1]=0
# C3


# In[12]:


# np.hstack((C, C))


# In[13]:


# for i in range(len(authors_index):
    # for j in 


# In[14]:


# authors_index


# # Loop over papers to create edgelist

# In[ ]:





# # EXTRA CODE

# In[ ]:





# In[249]:


type(authors_index)


# In[251]:


len(authors_index)


# In[241]:


mylist = authors_index[:10]


# In[242]:


mylist


# In[244]:


mylist.sort()


# In[246]:


# sort(mylist)
cars = ['Ford', 'BMW', 'Volvo']

cars.sort(reverse=True)


# In[247]:


cars


# In[ ]:




