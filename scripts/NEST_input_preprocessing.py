import numpy as np
import csv
import pickle
from scipy import sparse
import scipy.io as sio
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import stlearn as st
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
from typing import List
import qnorm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict
import pandas as pd
import gzip
from kneed import KneeLocator
import copy 
#sys.path.append("/home/gw/code/utility/altairThemes/")
#if True:  # In order to bypass isort when saving
#    import altairThemes
import altairThemes
import altair as alt
import argparse

spot_diameter = 89.43 #pixels


####################### Set the name of the sample you want to visualize ###################################

data_name = 'PDAC_64630' #LUAD_GSM5702473_TD1


if data_name == 'PDAC_64630':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' , help='The path to dataset') 
    parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset')
    #parser.add_argument( '--model_name', type=str, default='gat_2attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()


####### get the gene id, cell barcode, cell coordinates ######

if data_name == 'PDAC_64630':
    adata_h5 = st.Read10X(path=args.data_path, count_file='filtered_feature_bc_matrix.h5') #count_file=args.data_name+'_filtered_feature_bc_matrix.h5' )
    print(adata_h5)
    sc.pp.filter_genes(adata_h5, min_cells=1)
    print(adata_h5)
    gene_ids = list(adata_h5.var_names)
    coordinates = adata_h5.obsm['spatial']
    cell_barcode = np.array(adata_h5.obs.index)
    temp = adata_h5.X


######################### get the cell vs gene matrix ##################
'''
temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(temp)))  
adata_X = np.transpose(temp)  
cell_vs_gene = copy.deepcopy(adata_X)
'''
##################### make cell metadata: barcode_info ###################################
 
i=0
barcode_serial = dict()
for cell_code in cell_barcode:
    barcode_serial[cell_code]=i
    i=i+1
    
i=0
barcode_info=[]
for cell_code in cell_barcode:
    barcode_info.append([cell_code, coordinates[i,0],coordinates[i,1], 0]) # last entry will hold the component number later
    i=i+1

'''
i=0
node_id_sorted_xy=[]
for cell_code in cell_barcode:
    node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
    i=i+1
	
node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
'''


####### load annotations ##############################################

if data_name == 'PDAC_64630':
    pathologist_label_file='/cluster/home/t116508uhn/IX_annotation_artifacts.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)	
    	
    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]

