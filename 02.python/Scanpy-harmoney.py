#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc


# In[2]:


import random
random.seed(0)
import numpy as np 
np.random.seed(0)


# In[3]:


import scanpy.external as sce
import harmonypy


# In[4]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[5]:


results_file = './human_bone_marrow.h5ad'  # the file that will store the analysis results


# In[6]:


import os
path = '../02.downladdata_human-2022.7.20/08.result_file/'
file_dir = os.listdir(path)
adlist=[]
for file in file_dir:
    basename="human_batch1."+str.split(file,'.')[0]
    #print(basename)
    pathfile=path+file
    single_adata=sc.read_h5ad(pathfile)
    single_adata.var_names_make_unique()
    single_adata.obs_names_make_unique()
    single_adata.obs['batch']=basename
    adlist.append(single_adata)


# In[7]:


import os
path = '../04.human_download-2022.8.27/03.dynamo/result/'
file_dir = os.listdir(path)
for file in file_dir:
    basename="human_batch2."+file
    pathfile=path+file
    single_adata=sc.read_h5ad(pathfile)
    single_adata.var_names_make_unique()
    single_adata.obs_names_make_unique()
    # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    single_adata.obs['batch']=basename
    adlist.append(single_adata)


# In[8]:


import os
path = '../06.human-download-2022.9.8/04.picture/result/'
file_dir = os.listdir(path)
for file in file_dir:
    basename="human_batch3."+file
    pathfile=path+file
    single_adata=sc.read_h5ad(pathfile)
    single_adata.var_names_make_unique()
    single_adata.obs_names_make_unique()
    # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    single_adata.obs['batch']=basename
    adlist.append(single_adata)


# In[9]:


adlist


# In[10]:


import anndata as ad
adata = ad.concat(adlist,join='outer')


# In[11]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`


# In[12]:


adata.obs_names_make_unique()


# In[13]:


adata


# In[14]:


adata.obs


# In[15]:


adata.X


# # Preprocessing

# In[16]:


sc.pl.highest_expr_genes(adata, n_top=20, )


# In[17]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[18]:


adata.var['mt'] = adata.var_names.str.startswith('MT')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[19]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[20]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[21]:


#adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#adata = adata[adata.obs.pct_counts_mt < 5, :]


# In[22]:


sc.pp.normalize_total(adata, target_sum=1e4)


# In[23]:


sc.pp.log1p(adata)


# In[24]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[25]:


sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)


# # Principal component analysis

# In[26]:


sc.tl.pca(adata,svd_solver='arpack',random_state=0)


# In[27]:


sc.pp.neighbors(adata, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.1)


# In[28]:


sc.pl.scatter(adata, basis='umap', color='batch')


# In[29]:


adata


# In[30]:


sc.pl.scatter(adata, basis='umap', color='majority_voting')


# In[31]:


sce.pp.harmony_integrate(adata,'batch',theta=2,lamb=0.5)


# In[32]:


adata.var


# In[33]:


adata.write(results_file)


# In[34]:


adata


# # Computing the neighborhood graph

# In[35]:


sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40,random_state=0,use_rep="X_pca_harmony")


# # Embedding the neighborhood graph

# In[36]:


sc.tl.umap(adata,random_state=0, min_dist=0.2
          )


# In[37]:


sc.pl.umap(adata, color=['batch', 'GATA1', 'CST3','CD3E'])


# In[39]:


sc.pl.scatter(adata, basis='umap', color='majority_voting')


# In[40]:


HSC_markers=['CD34','KIT','IL7R', 'ATXN1', 'FLT3']


# In[41]:


sc.pl.scatter(adata, basis='umap', color=HSC_markers)


# # Clustering the neighborhood graph

# In[42]:


sc.tl.leiden(adata,random_state=0)


# In[43]:


get_ipython().run_line_magic('pinfo', 'sc.pl.umap')


# In[44]:


sc.pl.umap(adata, color=['leiden'])


# In[45]:


sc.pl.umap(adata, color=['leiden'],legend_loc='on data')


# In[46]:


#adata.obs['class']=adata.obs['leiden']
tmp={}
#color_key={}
for i in range(0,(len(adata.obs))):
    if adata.obs['batch'].str.contains('human_batch1')[i]:
        tmp[i]= adata.obs['batch'][i]
    else:
        tmp[i]='mon'
adata.obs['batch_class']=list(tmp.values())


# In[47]:


adata.obs


# In[48]:


color_key={
    'human_batch1.SRR12471880':'#336666',
    'human_batch1.SRR12471881':'#CCFFFF',
    'human_batch1.SRR12471882':'#FFCCCC',
    'human_batch1.SRR12471883':'#CCCCFF',
    'human_batch1.SRR12471884':'#990033',
    'human_batch1.SRR12471885':'#CCFF99',
    'human_batch1.SRR12471886':'#99CCCC',
    'human_batch1.SRR12471887':'#003366',
    'human_batch1.SRR12471888':'#666699',
    'mon':'#CCCCCC'
}


# In[49]:


adata.obs['batch_class']


# In[50]:


get_ipython().run_line_magic('pinfo', 'sc.pl.umap')


# In[51]:


adata2=adata[adata.obs['batch_class'].str.contains('human_batch1')]


# In[52]:


sc.pl.umap(adata, color='batch_class',palette=color_key)


# In[53]:


sc.pl.umap(adata2, color='batch_class',palette=color_key)


# In[54]:


tmp_index=np.arange(41415,158555).tolist()+np.arange(1,41415).tolist()


# In[55]:


tmp_index


# In[56]:


adata3=adata[tmp_index]


# In[57]:


sc.pl.umap(adata3, color=['batch_class'],palette=color_key)


# In[58]:


len(adata[adata.obs['batch'].str.contains('human_batch2')])


# In[59]:


41415+33434


# In[60]:


#adata.obs['class']=adata.obs['leiden']
tmp={}
#color_key={}
for i in range(0,(len(adata.obs))):
    if adata.obs['batch'].str.contains('human_batch2')[i]:
        tmp[i]= adata.obs['batch'][i]
    else:
        tmp[i]='mon'
adata.obs['batch_class2']=list(tmp.values())
tmp_index=np.arange(41415,74849).tolist()+np.arange(1,41415).tolist()+np.arange(41415,158555).tolist()
adata4=adata[tmp_index]


# In[61]:


sc.pl.umap(adata4, color=['batch_class2'])


# In[ ]:


#adata.obs['class']=adata.obs['leiden']
tmp={}
#color_key={}
for i in range(0,(len(adata.obs))):
    if adata.obs['batch'].str.contains('human_batch3')[i]:
        tmp[i]= adata.obs['batch'][i]
    else:
        tmp[i]='mon'
adata.obs['batch_class3']=list(tmp.values())
tmp_index=np.arange(41415,158555).tolist()+np.arange(1,41415).tolist()+np.arange(41415,74849).tolist()
adata5=adata[tmp_index]


# In[ ]:


sc.pl.umap(adata5, color=['batch_class3'])


# In[ ]:


sc.pl.umap(adata, color=['MKI67'],legend_loc='on data')


# In[ ]:


adata


# In[ ]:


adata.obs


# In[ ]:


#adata.obs['class']=adata.obs['leiden']
tmp={}
for i in range(0,(len(adata.obs))):
   tmp[i]= adata.obs['leiden'][i]+":"+adata.obs['majority_voting'][i]

adata.obs['class']=list(tmp.values())


# In[ ]:


adata.obs


# In[ ]:


tmp1=adata[adata.obs['leiden']=='0']


# In[ ]:


tmp2=adata[adata.obs['batch'].str.contains('human_batch1')]


# In[ ]:


len(tmp2.obs['class'])


# In[ ]:


set(tmp1.obs['class'])


# In[ ]:


sc.pl.umap(adata, color=['majority_voting'])


# In[ ]:


adata.obs['leiden']


# In[ ]:


sc.pl.umap(adata, color=HSC_markers)


# In[ ]:


CLP_marker=['CD34','CD38']


# In[ ]:


pro-B marker=['CD34','CD38','TdT', 'CD179B' 'Vpreb' 'CD10' 'ILTRA']


# In[ ]:


sc.pl.umap(adata, color=CLP_marker)


# In[ ]:


sc.pl.umap(adata, color=pro-B marker)


# In[ ]:


adata.write(results_file)


# # Celltyplist

# In[56]:


import celltypist
from celltypist import models


# In[57]:


# Enabling `force_update = True` will overwrite existing (old) models.
models.download_models(force_update = True)
models.models_description()


# In[58]:


# Indeed, the `model` argument defaults to `Immune_All_Low.pkl`.
model = models.Model.load(model = 'Immune_All_Low.pkl')
model.cell_types


# In[59]:


predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)


# In[60]:


predictions.predicted_labels


# In[61]:


# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata_new = predictions.to_adata()


# In[62]:


get_ipython().run_line_magic('pinfo', 'sc.pl.umap')


# In[63]:


sc.pl.umap(adata, color = ['majority_voting'], legend_loc = 'right margin',legend_fontsize='12')


# In[64]:


adata_color=pd.read_table('../04.human_download-2022.8.27/04.h5ad/color_v3.list')
color_key={}
for i in range(0,len(model.cell_types)):
    color_key[model.cell_types[i]]=adata_color['number'][i]
color_key['Classical monocytes']='#1f77b4'
color_key['HSC/MPP']='#d62728'
color_key['MEMP']='#aa40fc'
color_key['Erythrophagocytic macrophages']='#ff7f0e'
color_key['Mast cells']='#8c564b'
color_key['Megakaryocyte precursor']='#e377c2'
color_key['Neutrophil-myeloid progenitor']='#17becf'
color_key['Large pre-B cells']='#98df8a'
color_key['CMP']='#ffbb78'
color_key['Pro-B cells']='#ff9896'
color_key['CD16-NK-cells']='#279e68'


# In[65]:


sc.pl.umap(adata, color = ['majority_voting'], legend_loc = 'right margin',legend_fontsize='12',palette=color_key)


# In[66]:


CLP_markers=['CD34','CD38','CD117', 'CD1798']


# In[67]:


HSC_markers=['CD34','KIT','IL7R', 'ATXN1', 'FLT3']


# In[68]:


sc.pl.umap(adata, color=HSC_markers)


# # dynamo

# In[69]:


import dynamo as dyn 
dyn.get_all_dependencies_version()
dyn.configuration.set_figure_params('dynamo', background='white')


# In[70]:


dyn.pp.recipe_monocle(adata)
dyn.tl.dynamics(adata, model='stochastic', cores=3)
dyn.tl.reduceDimension(adata)


# In[71]:


dyn.pl.umap(adata, color='majority_voting',alpha=1,pointsize=0.01)


# In[72]:


# dyn.tl.gene_wise_confidence(adata, group='group', lineage_dict={'Progenitor': ['terminal_cell_state']})
dyn.tl.gene_wise_confidence(adata, group='majority_voting', lineage_dict={'HSC/MPP': ['Early erythroid']})
dyn.pl.phase_portraits(adata, genes=adata.var_names[adata.var.use_for_dynamics][:4], figsize=(6, 4), color='majority_voting')
dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'})


# In[73]:


dyn.tl.cell_wise_confidence(adata)
# dyn.tl.confident_cell_velocities(adata, group='group', lineage_dict={'Progenitor': ['terminal_cell_state']},)
dyn.tl.confident_cell_velocities(adata, group='majority_voting', lineage_dict={'HSC/MPP': ['Early erythroid']})


# In[74]:


dyn.pl.streamline_plot(adata, color=['majority_voting'], color_key = color_key,pointsize=0.5, alpha=1,
                       basis='umap', show_legend='on data', show_arrowed_spines=True)


# # Finding marker genes

# In[75]:


sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[76]:


sc.tl.rank_genes_groups(adata, 'majority_voting', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[77]:


sc.settings.verbosity = 2  # reduce the verbosity


# In[83]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[84]:


sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[ ]:


marker_genes = ["CD34","KIT","IL7R","ATXN1","FLT3"]


# In[ ]:


adata = sc.read(results_file)


# In[ ]:


pd.DataFrame(adata.uns['rank_genes_groups']['names'])


# In[ ]:


result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})


# In[ ]:


sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)


# In[ ]:


adata = sc.read(results_file)


# In[ ]:


adata.obs


# In[ ]:


sc.pl.violin(adata, marker_genes, groupby='majority_voting')


# In[ ]:





# In[ ]:


sc.pl.violin(adata,marker_genes, groupby='louvain')


# In[ ]:


adata.write(results_file)


# In[ ]:





# In[ ]:




