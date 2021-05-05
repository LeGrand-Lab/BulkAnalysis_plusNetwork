# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import holoviews, pygraphviz
import networkx as nx
import pandas as pd
import matplotlib.patches as mpatches

myroot = "~/bulk_analysis/"

print(os.getcwd())
natmiVizOut = "natmiD7/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/"
tb = pd.read_csv(f'{myroot}{natmiVizOut}Edges.csv',sep=",", header=0)

toy = tb.sample(n=10)



# make UNIQUE the ligand and receptor names by combining celltype alongside
toy['uniq_Ligand_symbol'] = toy['Ligand symbol'] + toy['Sending cluster']
toy['uniq_Receptor_symbol'] = toy['Receptor symbol'] + toy['Target cluster']

G = nx.DiGraph(age_day="Young_D7")
G.add_node('mooA', celltype = 'A', genesym = 'uuu', color='red',specif=0.01)
G.add_node('mooB', celltype= 'B', genesym = 'uuu', color='blue',specif=0.2)
G.add_edge('mooA','mooB', origtype='A')
G.add_edge('mooA','mooA', origtype='A')
node_labels = nx.get_node_attributes(G,'omg')


pos = nx.circular_layout(G)
nodes = G.nodes()
ccc = [nodes[n]['color'] for n in nodes]
nx.draw(G,pos)
selfloops = [('mooA','mooA')]
nx.draw_networkx_labels(G, pos, labels=node_labels)
nx.draw_networkx_nodes(G, pos, node_color=ccc)
nx.draw_networkx_edges(G, pos, edgelist=selfloops, arrowstyle="<|-", style="dashed")
plt.savefig('pretoy.png')
plt.show()

# =================
from networkx.drawing.nx_agraph import write_dot 
write_dot(G,'graph.dot')
# =================

labels = ['Label 1', 'Label 2']
colors = ['cyan','gray']
 
fig = plt.figure(figsize=(2, 1.25))
patches = [
    mpatches.Patch(color=color, label=label)
    for label, color in zip(labels, colors)]
fig.legend(patches, labels, loc='center', frameon=False)
plt.show()


# !! ??
import matplotlib.pyplot as plt
plt.ion()
import networkx
import netgraph # pip install netgraph

# Construct sparse, directed, weighted graph
total_nodes = 20
weights = np.random.rand(total_nodes, total_nodes)
connection_probability = 0.1
is_connected = np.random.rand(total_nodes, total_nodes) <= connection_probability
graph = np.zeros((total_nodes, total_nodes))
graph[is_connected] = weights[is_connected]

# construct a networkx graph
g = networkx.from_numpy_array(graph, networkx.DiGraph)

# decide on a layout
pos = networkx.layout.spring_layout(g)

# Create an interactive plot.
# NOTE: you must retain a reference to the object instance!
# Otherwise the whole thing will be garbage collected after the initial draw
# and you won't be able to move the plot elements around.
plot_instance = netgraph.InteractiveGraph(graph, node_positions=pos)

######## drag nodes around #########

# To access the new node positions:
node_positions = plot_instance.node_positions