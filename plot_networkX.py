import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse 
from matplotlib import rc, cm
import networkx as nx
# import biographs as bg
import os
import plotly.graph_objects as go
# from graphein.protein.config import ProteinGraphConfig
# from graphein.protein.visualisation import plotly_protein_structure_graph
# from graphein.protein.graphs import construct_graph
# from graphein.protein.analysis import plot_edge_type_distribution, graph_summary
# from graphein.protein.edges.distance import (
#     add_aromatic_interactions,
#     add_disulfide_interactions,
#     add_hydrophobic_interactions,
#     add_peptide_bonds,
# )

# sns.set_context("notebook")
# plt.style.use("dark_background")
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='CORRELATION FILE')
parser.add_argument('--pdb', type=str, default='',
                    help='PDB FILE with corrected CHAIN')
args = parser.parse_args()

ifile = args.input
pdb_file = args.pdb


with open(ifile ,'r') as f:
    lines = f.readline()    # SKIP THE FIRST LINE
    # print(lines)
    meta_info = lines.split(' ')

    rows = int(meta_info[0])
    cols = int(meta_info[2])
    data = []
    lines = f.readlines()
    for line in lines:
        data.extend(line.split())

    data.pop(-1)    # REMOVE THE LAST ']'
    # print(len(data))
    data = np.array(data,dtype=np.float64)
    matrix = np.reshape(data,(rows,cols))

resid_a = np.arange(73,928)
resid_b = np.arange(73,928)
resid_list = np.concatenate((resid_a,resid_b))
df = pd.DataFrame(matrix,columns=resid_list,index=resid_list)
df = df.fillna(0)
# print(df)

G=nx.Graph(name='Protein Interaction Graph')
interactions = np.array(df)
for i in range(len(interactions)):
    interaction = interactions[i]
    # print(interaction)
    a = interaction[0] # protein a node
    b = interaction[1] # protein b node
    w = float(interaction[2]) # score as weighted edge where high scores = low weight
    G.add_weighted_edges_from([(a,b,w)]) # add weighted edge to graph

# pos = nx.spring_layout(G) # position the nodes using the spring layout
# plt.figure(figsize=(11,11),facecolor=[0.7,0.7,0.7,0.4])
# nx.draw_networkx(G)
# plt.axis('off')
# plt.show()

# function to rescale list of values to range [newmin,newmax]
def rescale(l,newmin,newmax):
    arr = list(l)
    return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]
# use the matplotlib plasma colormap
graph_colormap = cm.get_cmap('plasma', 12)
# node color varies with Degree
c = rescale([G.degree(v) for v in G],0.0,0.9) 
c = [graph_colormap(i) for i in c]
# node size varies with betweeness centrality - map to range [10,100] 
bc = nx.betweenness_centrality(G) # betweeness centrality
s =  rescale([v for v in bc.values()],1500,7000)
# edge width shows 1-weight to convert cost back to strength of interaction 
ew = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,4)
# edge color also shows weight
ec = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,1)
ec = [graph_colormap(i) for i in ec]
# view raw

pos = nx.spring_layout(G)
plt.figure(figsize=(19,9),facecolor=[0.7,0.7,0.7,0.4])
nx.draw_networkx(G, pos=pos, with_labels=True, node_color=c, node_size=s,edge_color= ec,width=ew,
                 font_color='white',font_weight='bold',font_size='9')
plt.axis('off')
plt.show()