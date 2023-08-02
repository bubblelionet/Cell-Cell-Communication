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
##########################################################
# readCsv, preprocessDf, plot: these three functions are taken from GW's repository                                                                                                                                                                     /mnt/data0/gw/research/notta_pancreatic_cancer_visium/plots/fatema_signaling/hist.py                                                                                                                                                                                         
import scipy.stats


def readCsv(x):
  """Parse file."""
  #colNames = ["method", "benchmark", "start", "end", "time", "memory"]
  df = pd.read_csv(x, sep=",")

  return df

def preprocessDf(df):
  """Transform ligand and receptor columns."""
  df["ligand-receptor"] = df["ligand"] + '-' + df["receptor"]
  df["component"] = df["component"] #.astype(str).str.zfill(2)

  return df


def plot(df):
  set1 = altairThemes.get_colour_scheme("Set1", len(df["component"].unique()))
  set1[0] = '#000000'
  base = alt.Chart(df).mark_bar().encode(
            x=alt.X("ligand-receptor:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y=alt.Y("count()"),
            color=alt.Color("component:N", scale = alt.Scale(range=set1)),
            order=alt.Order("component:N", sort="ascending"),
            tooltip=["component"]
        )
  p = base

  return p

####################### Set the name of the sample you want to visualize ###################################

data_name = 'PDAC_64630' #LUAD_GSM5702473_TD1



##########################################################
if data_name == 'LUAD_GSM5702473_TD1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='LUAD_GSM5702473_TD1', help='The name of dataset')
    parser.add_argument( '--model_name', type=str, default='gat_r1_3attr', help='model name')
    
    args = parser.parse_args()  	
#############################################################   
elif data_name == 'V1_Human_Lymph_Node_spatial':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='V1_Human_Lymph_Node_spatial', help='The name of dataset')
    parser.add_argument( '--model_name', type=str, default='gat_r1_2attr', help='model name')
    parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()


elif data_name == 'PDAC_64630':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset')
    #parser.add_argument( '--model_name', type=str, default='gat_2attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()


elif data_name == 'PDAC_140694':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/V10M25-60_C1_PDA_140694_Pa_P_Spatial10x/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_140694', help='The name of dataset')
    args = parser.parse_args()

####### get the gene id, cell barcode, cell coordinates ######

if data_name == 'LUAD_GSM5702473_TD1':
    # read the mtx file
    temp = sc.read_10x_mtx(args.data_path)
    print(temp)
    sc.pp.filter_genes(temp, min_cells=1)
    print(temp)
    
    gene_ids = list(temp.var_names) 
    cell_barcode = np.array(temp.obs.index)
    
    # now read the tissue position file. It has the format: 
    #df = pd.read_csv('/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/spatial/tissue_positions_list.csv', sep=",",header=None)   # read dummy .tsv file into memory
    df = pd.read_csv('/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv', sep=",",header=None)   # read dummy .tsv file into memory
    tissue_position = df.values
    barcode_vs_xy = dict() # record the x and y coord for each spot
    for i in range (0, tissue_position.shape[0]):
        barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][5], tissue_position[i][4]] #for some weird reason, in the .h5 format, the x and y are swapped
        #barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][4], tissue_position[i][5]] 
    
    coordinates = np.zeros((cell_barcode.shape[0], 2)) # insert the coordinates in the order of cell_barcodes
    for i in range (0, cell_barcode.shape[0]):
        coordinates[i,0] = barcode_vs_xy[cell_barcode[i]][0]
        coordinates[i,1] = barcode_vs_xy[cell_barcode[i]][1]
    

else:
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
if data_name == 'LUAD_GSM5702473_TD1':
    '''
    ccc_too_many_cells_LUAD = pd.read_csv('/cluster/projects/schwartzgroup/fatema/CCST/exp2_D1_ccc_toomanycells_cluster.csv')
    ccc_too_many_cells_LUAD_dict = dict()
    for i in range(0, len(ccc_too_many_cells_LUAD)):
        ccc_too_many_cells_LUAD_dict[ccc_too_many_cells_LUAD['cell'][i]] = int(ccc_too_many_cells_LUAD['cluster'][i])
    
    for i in range(0, len(barcode_info)):
        barcode_info[i][3] = ccc_too_many_cells_LUAD_dict[barcode_info[i][0]]
    	
    '''
    barcode_type=dict()
    for i in range (0, len(barcode_info)):
        barcode_type[barcode_info[i][0]] = 0 
    
    
elif data_name == 'V1_Human_Lymph_Node_spatial':
    '''
    pathologist_label_file='/cluster/home/t116508uhn/human_lymphnode_Spatial10X_manual_annotations.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        if pathologist_label[i][1] == 'GC': 
            barcode_type[pathologist_label[i][0]] = 1 
        #elif pathologist_label[i][1] =='':
        #    barcode_type[pathologist_label[i][0]] = 0 #'stroma_deserted'
        else:
            barcode_type[pathologist_label[i][0]] = 0
    '''
    pathologist_label_file='/cluster/home/t116508uhn/64630/spot_vs_type_dataframe_V1_HumanLympNode.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    	
    spot_label = []
    for i in range (1, len(pathologist_label)):
        spot_label.append([pathologist_label[i][0], float(pathologist_label[i][1]), float(pathologist_label[i][2]), float(pathologist_label[i][3])])
        
    spot_label = sorted(spot_label, key = lambda x: x[3], reverse=True) # descending order of 
    barcode_Tcell = []
    barcode_B = []
    barcode_GC = []
    for i in range (0, len(spot_label)):
        if spot_label[i][1] >= (spot_label[i][2] + spot_label[i][3])*2:
            barcode_Tcell.append(spot_label[i][0])
        elif spot_label[i][2] >= (spot_label[i][1] + spot_label[i][3])*2:
            barcode_B.append(spot_label[i][0])
        elif spot_label[i][3] >= (spot_label[i][1] + spot_label[i][2])*2:
            barcode_GC.append(spot_label[i][0])
           
    barcode_type=dict()
    for i in range (0, len(barcode_Tcell)):
        barcode_type[barcode_Tcell[i]] = 1 # tcell
        
    for i in range (0, len(barcode_B)):
        barcode_type[barcode_B[i]] = 2
      
    for i in range (0, len(barcode_GC)):
        barcode_type[barcode_GC[i]] = 3
        
    for i in range (0, len(spot_label)):
        if spot_label[i][0] not in barcode_type:
            barcode_type[spot_label[i][0]] = 0
    #############################################################################
    '''
    data_list=dict()
    data_list['pathology_label']=[]
    data_list['X']=[]
    data_list['Y']=[]
    
    for i in range (0, len(barcode_info)):
        data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(barcode_info[i][2])
    
    
    data_list_pd = pd.DataFrame(data_list)
    set1 = altairThemes.get_colour_scheme("Set1", 4)
    set1[0] = '#000000'
    
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
        color=alt.Color('pathology_label:N', scale=alt.Scale(range=set1)),
        tooltip=['pathology_label']
    )
    
    save_path = '/cluster/home/t116508uhn/64630/'
    chart.save(save_path+'V1_humanLymphNode.html') #   
    '''            
elif data_name == 'PDAC_64630':
    pathologist_label_file='/cluster/home/t116508uhn/IX_annotation_artifacts.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)	
    	
    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]
        
    '''
    ########## sabrina ###########################################	
    pathologist_label_file='/cluster/projects/schwartzgroup/fatema/find_ccc/singleR_spot_annotation_Sabrina.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        if pathologist_label[i][1] == 'tumor': #'Tumour':
            barcode_type[pathologist_label[i][0]] = 1 #'tumor'
        elif pathologist_label[i][1] =='stroma_deserted':
            barcode_type[pathologist_label[i][0]] = 0 #'stroma_deserted'
        elif pathologist_label[i][1] =='acinar_reactive':
            barcode_type[pathologist_label[i][0]] = 2 #'acinar_reactive'
        else:
            barcode_type[pathologist_label[i][0]] = 0 #'zero' 
    '''	
    #################################################################		
elif data_name == 'PDAC_140694':
    spot_type = []
    pathologist_label_file='/cluster/home/t116508uhn/V10M25-060_C1_T_140694_Histology_annotation_IX.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
            spot_type.append(line[1])
            
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        if 'tumor_LVI' in pathologist_label[i][1]:
            barcode_type[pathologist_label[i][0]] = 'tumor_LVI'
        elif 'tumor_PNI' in pathologist_label[i][1]:
            barcode_type[pathologist_label[i][0]] = 'tumor_PNI'
        elif 'tumor_stroma' in pathologist_label[i][1]:
            barcode_type[pathologist_label[i][0]] = 'tumor_stroma'
        elif 'tumor_vs_acinar' in pathologist_label[i][1]:
            barcode_type[pathologist_label[i][0]] = 'tumor_vs_acinar'     
        elif 'nerve' in pathologist_label[i][1]: 
            barcode_type[pathologist_label[i][0]] = 'nerve'
        elif 'vessel' in pathologist_label[i][1]: 
            barcode_type[pathologist_label[i][0]] = 'vessel'
        #elif pathologist_label[i][1] =='stroma_deserted':
        #    barcode_type[pathologist_label[i][0]] = 'stroma_deserted'
        #elif pathologist_label[i][1] =='acinar_reactive':
        #    barcode_type[pathologist_label[i][0]] = 'acinar_reactive'
        else:
            barcode_type[pathologist_label[i][0]] = 'others'
           


######################### read the NEST output in csv format ####################################################


filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
inFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]
df = pd.read_csv(inFile, sep=",")
csv_record_final = df.values.tolist()
df_column_names = list(df.columns)
csv_record_final = [df_column_names] + csv_record_final

# columns are: from_cell, to_cell, ligand_gene, receptor_gene, attention_score, component, from_id, to_id

################################################################################

## change the csv_record_final here if you want histogram for specific components/regions only. e.g., if you want to plot only stroma region, or tumor-stroma regions etc.    ##
'''
region_of_interest = [...] 
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    # if both of ligand and receptors are tumors, or both of them are non-tumors, then remove it. Because we want to see what ccc is happening between tumor and non-tumor. 
    if ((barcode_type[csv_record_final[record_idx][0]] == 'tumor' and barcode_type[csv_record_final[record_idx][1]] == 'tumor') or (barcode_type[csv_record_final[record_idx][0]] != 'tumor' and barcode_type[csv_record_final[record_idx][1]] != 'tumor')):
        csv_record_final[record_idx][5] = 0 # label it 0 so that it is not considered during ploting and making histogram
	
'''

###########	dictionary of those spots who are participating in CCC ##################
active_spot = defaultdict(list)
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    record = csv_record_final[record_idx]
    i = record[6]
    pathology_label = barcode_type[barcode_info[i][0]]
    component_label = record[5]
    X = barcode_info[i][1]
    Y = -barcode_info[i][2]
    opacity = record[4]
    active_spot[i].append([pathology_label, component_label, X, Y, opacity])
    
    j = record[7]
    pathology_label = barcode_type[barcode_info[j][0]]
    component_label = record[5]
    X = barcode_info[j][1]
    Y = -barcode_info[j][2]
    opacity = record[4]   
    active_spot[j].append([pathology_label, component_label, X, Y, opacity])
    ''''''
    
######### color the spots in the plot with opacity = attention score #################
opacity_list = []
for i in active_spot:
    sum_opacity = []
    for edges in active_spot[i]:
        sum_opacity.append(edges[4])
        
    avg_opacity = np.max(sum_opacity) #np.mean(sum_opacity)
    opacity_list.append(avg_opacity) 
    active_spot[i]=[active_spot[i][0][0], active_spot[i][0][1], active_spot[i][0][2], active_spot[i][0][3], avg_opacity]

min_opacity = np.min(opacity_list)
max_opacity = np.max(opacity_list)
min_opacity = min_opacity - 5

################## altair plot #####################################################################
data_list=dict()
data_list['pathology_label']=[]
data_list['component_label']=[]
data_list['X']=[]
data_list['Y']=[]   
data_list['opacity']=[] 

for i in range (0, len(barcode_info)):        
    if i in active_spot:
        data_list['pathology_label'].append(active_spot[i][0])
        data_list['component_label'].append(active_spot[i][1])
        data_list['X'].append(active_spot[i][2])
        data_list['Y'].append(active_spot[i][3])
        data_list['opacity'].append((active_spot[i][4]-min_opacity)/(max_opacity-min_opacity))
        
    else:
        data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
        data_list['component_label'].append(0)
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(-barcode_info[i][2])
        data_list['opacity'].append(0.1)


#id_label= len(list(set(data_list['component_label'])))
data_list_pd = pd.DataFrame(data_list)
set1 = altairThemes.get_colour_scheme("Set1", id_label)
set1[0] = '#000000'
chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
    alt.X('X', scale=alt.Scale(zero=False)),
    alt.Y('Y', scale=alt.Scale(zero=False)),
    shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
    color=alt.Color('component_label:N', scale=alt.Scale(range=set1)),
    #opacity=alt.Opacity('opacity:N'), #"opacity",
    tooltip=['component_label'] #,'opacity'
)#.configure_legend(labelFontSize=6, symbolLimit=50)

# output 6
save_path = '/cluster/home/t116508uhn/64630/'
chart.save(save_path+'altair_plot_test.html')
###################################  Histogram plotting #################################################################################
filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")
inFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]



df = readCsv(inFile)
df = preprocessDf(df)
outPathRoot = inFile.split('.')[0]
p = plot(df)
outPath = '/cluster/home/t116508uhn/64630/test_hist_temp.html'
p.save(outPath)	# output 5
#####################################################################################################################




############################  Network Plot ######################
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
set1 = altairThemes.get_colour_scheme("Set1", id_label)
colors = set1
colors[0] = '#000000'
ids = []
x_index=[]
y_index=[]
colors_point = []
for i in range (0, len(barcode_info)):    
    ids.append(i)
    x_index.append(barcode_info[i][1])
    y_index.append(barcode_info[i][2])    
    colors_point.append(colors[barcode_info[i][3]]) 
  
max_x = np.max(x_index)
max_y = np.max(y_index)


from pyvis.network import Network
import networkx as nx

barcode_type=dict()
for i in range (1, len(pathologist_label)):
    if 'tumor'in pathologist_label[i][1]: #'Tumour':
        barcode_type[pathologist_label[i][0]] = 1
    else:
        barcode_type[pathologist_label[i][0]] = 0
    '''
    elif pathologist_label[i][1] == 'stroma_deserted':
        barcode_type[pathologist_label[i][0]] = 0
    elif pathologist_label[i][1] =='acinar_reactive':
        barcode_type[pathologist_label[i][0]] = 2
    else:
        barcode_type[pathologist_label[i][0]] = 'zero' #0
    '''
g = nx.MultiDiGraph(directed=True) #nx.Graph()
for i in range (0, len(barcode_info)):
    label_str =  str(i)+'_c:'+str(barcode_info[i][3])+'_'
    #if barcode_type[barcode_info[i][0]] == 'zero':
    #    continue
    if barcode_type[barcode_info[i][0]] == 0: #stroma
        marker_size = 'circle'
        label_str = label_str + 'stroma'
    elif barcode_type[barcode_info[i][0]] == 1: #tumor
        marker_size = 'box'
        label_str = label_str + 'tumor'
    else:
        marker_size = 'ellipse'
        label_str = label_str + 'acinar_reactive'
	
    g.add_node(int(ids[i]), x=int(x_index[i]), y=int(y_index[i]), label = label_str, pos = str(x_index[i])+","+str(-y_index[i])+" !", physics=False, shape = marker_size, color=matplotlib.colors.rgb2hex(colors_point[i]))    

nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True

count_edges = 0
for k in range (1, len(csv_record_final)):
    i = csv_record_final[k][6]
    j = csv_record_final[k][7]    
    ligand = csv_record_final[k][2]
    receptor = csv_record_final[k][3]
    title_str =  "L:"+ligand+", R:"+receptor
    edge_score = csv_record_final[k][4]
    g.add_edge(int(i), int(j), label = title_str, value=np.float64(edge_score), color=colors_point[i] ) 
    count_edges = count_edges + 1
     
nt.from_nx(g)
nt.show('mygraph.html')
cp mygraph.html /cluster/home/t116508uhn/64630/mygraph.html

from networkx.drawing.nx_agraph import write_dot
write_dot(g, "/cluster/home/t116508uhn/64630/test_interactive.dot")

