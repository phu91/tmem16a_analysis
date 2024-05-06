import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse 
from matplotlib import rc
import networkx as nx
# import biographs as bg
import os, re
import plotly.graph_objects as go
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.visualisation import plotly_protein_structure_graph
from graphein.protein.graphs import construct_graph
from graphein.protein.analysis import plot_edge_type_distribution, graph_summary
from graphein.protein.edges.distance import (
    add_aromatic_interactions,
    add_disulfide_interactions,
    add_hydrophobic_interactions,
    add_peptide_bonds,
)

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


def add_GC_edges(G: nx.Graph) -> nx.Graph:
    residues_node = [n for n, d in G.nodes(data=True)]
    for i, x in enumerate(residues_node):
        for j, y in enumerate(residues_node):
            gc = df_corr.iloc[i,j]
            G.add_edge(x, y, mutual_information=gc)
    return G
#### CONSTRUCT GENERALIZED CORRELATION MATRIX
with open(ifile ,'r') as f:
    lines = f.readline()    # SKIP THE FIRST LINE
    # print(lines)
    meta_info = lines.split(' ')

    rows = int(meta_info[0])
    cols = int(meta_info[2])
    long_matrix = []
    lines = f.readlines()
    for line in lines:
        long_matrix.extend(line.split())
    long_matrix.pop(-1)    # REMOVE THE LAST ']'
    # print(len(data))
    long_matrix = np.array(long_matrix,dtype=np.float64)
    # long_matrix[long_matrix <0.5 ] = 0
    matrix = np.reshape(long_matrix,(rows,cols))

resid_a = np.arange(73,928)
resid_b = np.arange(73,928)
resid_list = np.concatenate((resid_a,resid_b))
df_corr = pd.DataFrame(matrix,columns=resid_list,index=resid_list)
# df_chainA = df_corr.iloc[:855,:855] ## CHAIN A (before column 855)d
# df_chainB = df_corr.iloc[855:,855:] ## CHAIN B (after column 855)
# print(df_chainA.iloc[0,0])

### CONSTRUCT GRAPH

# Defined default configuration for graph
config = ProteinGraphConfig(granularity='CA')
# Add new edge function 
new_edge_funcs = {"edge_construction_functions": [add_GC_edges,add_aromatic_interactions]}
config = ProteinGraphConfig(**new_edge_funcs)

## Construct the protein (residue) graph
g = construct_graph(config=config,path=pdb_file,chain_selection=['A','B'])

# Filter edges based on 'kind' attribute (assuming 'kind' is an edge attribute)
filtered_edges = []
for u, v, a in g.edges(data=True):
    if a.get("mutual_information")<1000 and u!=v and a.get("distance")>11:
        filtered_edges.append((u, v))
        # print(u,v,a)

# ## CHECKING 
# for u, v, a in g.edges(data=True):
#     print(u,v,a)

# # Create a subgraph with only the filtered edges
filtered_subgraph = g.edge_subgraph(filtered_edges)
# filtered_subgraph = g
# u = graph_summary(filtered_subgraph)
# u.to_csv("NODE_features.csv")

###Plot the graph with only the specified 'kind' of edges
p = plotly_protein_structure_graph(
    filtered_subgraph,
    colour_edges_by="mutual_information",
    colour_nodes_by="chain",
    label_node_ids=False,
    # node_alpha=0.5,
    # node_size_min=10,
    # node_size_feature="betweenness_centrality",
    plot_title="Protein graph created using a user-defined function that connects all proline.",
    # figsize=(1792, 1120)
)
p.show()
### MISCELLANEOUS ###
# plt.locator_params(axis='y', nbins=30)
# plt.locator_params(axis='x', nbins=30)
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# # plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
# plt.show()
