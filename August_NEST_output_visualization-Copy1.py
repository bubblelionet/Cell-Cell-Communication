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



######################### read the NEST output in csv format ####################################################


filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
inputFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]
df = pd.read_csv(inputFile, sep=",")
csv_record_final = df.values.tolist()
df_column_names = list(df.columns)
csv_record_final = [df_column_names] + csv_record_final

# columns are: from_cell, to_cell, ligand_gene, receptor_gene, attention_score, component, from_id, to_id

################change csv_record_final for specific components/regions only ##############
# by specifying regions we DO NOT WANT TO PLOT ####

# i.e. plot only stroma region/ tumor-stroma regions. 

region_of_interest = [...] 
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, ignore.
    if ((barcode_type[csv_record_final[record_idx][0]] == 'tumor' and barcode_type[csv_record_final[record_idx][1]]== 'tumor') or (barcode_type[csv_record_final[record_idx][0]] != 'tumor' and barcode_type[csv_record_final[record_idx][1]] != 'tumor')):csv_record_final[record_idx][5]= 0 # This will exclude all tumor-tumor and stroma-stroma interactions
        
# i.e. plot only for CCC between one component and everything else, i.e. set ligand as component 13 and receptor as everything else but 13
# i.e. plot only within a certain region by setting a boundary xy coordinate



########### create active_spot file: dictionary of spots participating in CCC ##################
active_spot = defaultdict(list)

for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    record = csv_record_final[record_idx]
    
    # ligand
    i = record[6] # correspond to the column number in the input file 
    pathology_label = barcode_type[barcode_info[i][0]] 
    component_label = record[5]
    X = barcode_info[i][1] # barcode_info file was created earlier
    Y = -barcode_info[i][2]
    opacity = record[4]
    active_spot[i].append([pathology_label, component_label, X, Y, opacity])
    
    # receptor
    j = record[7] 
    pathology_label = barcode_type[barcode_info[j][0]]
    component_label = record[5]
    X = barcode_info[j][1]
    Y = -barcode_info[j][2]
    opacity = record[4]   
    active_spot[j].append([pathology_label, component_label, X, Y, opacity])
    
    
######### color the spots in the plot with opacity = attention score #################
opacity_list = []
for i in active_spot:
    sum_opacity = []
    for edges in active_spot[i]:
        sum_opacity.append(edges[4])
        
    # avg_opacity = np.mean(sum_opacity)
    max_opacity = np.max(sum_opacity)
    opacity_list.append(max_opacity) 
    active_spot[i]=[active_spot[i][0][0], active_spot[i][0][1], active_spot[i][0][2], active_spot[i][0][3], max_opacity]

min_opacity = np.min(opacity_list)
max_opacity = np.max(opacity_list)
min_opacity = min_opacity - 5

################## altair plot #####################################################################

# put everything into pandas df
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

###########################  Histogram plotting ####################################################
filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")
inputFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]



df = readCsv(inputFile)
df = preprocessDf(df)
outPathRoot = inputFile.split('.')[0]
p = plot(df)
outPath = '/cluster/home/t116508uhn/64630/test_hist_temp.html'
p.save(outPath)	# output 5

#########################  Network Plot ######################
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


# visualize the actual communication
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

# need to convert dot format to pdf

