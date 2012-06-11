import networkx as nx
import matplotlib.pyplot as plt
import collections

def draw_coms(g, pos, cs):
    NODE_COLORS = 'bgcmy'
    NODE_SIZE = 500
    nx.draw_networkx_nodes(g, pos=pos,
                           node_color='0.5',
                           node_size=NODE_SIZE)
    vis = set()
    for i, c in enumerate(cs):
        dn = vis & c
        vis = vis | c
        nx.draw_networkx_nodes(g, pos=pos,
                               nodelist=c,
                               node_color=NODE_COLORS[i % len(NODE_COLORS)],
                               node_size=NODE_SIZE)
        nx.draw_networkx_nodes(g, pos=pos,
                               nodelist=dn,
                               node_color='r',
                               node_size=NODE_SIZE)
        nx.draw_networkx_edges(g, pos=pos,
                               width=2,
                               edge_color='0.2')
        nx.draw_networkx_labels(g, pos=pos, font_color='0.95', font_size=11)
        plt.show()

def draw_cores(g, pos, cn):
    NODE_COLORS = 'bgcmy'
    NODE_SIZE = 500
    nx.draw_networkx_nodes(g, pos=pos,
                           node_color=[NODE_COLORS[v-1] for v in cn.values()],
                           node_size=NODE_SIZE)
    nx.draw_networkx_edges(g, pos=pos,
                           width=2,
                           edge_color='0.2')
    nx.draw_networkx_labels(g, pos=pos, font_color='0.95', font_size=11)
    plt.show()
    

g = nx.read_gml('karate.gml')
#g = nx.read_edgelist('ca-AstroPh.txt')
cs = nx.k_clique_communities(g, 5)
pos = nx.spring_layout(g, iterations=200)
draw_coms(g, pos, cs)
cn = nx.core_number(g)
print cn
draw_cores(g, pos, cn)

#g = nx.read_edgelist('ca-AstroPh.txt', nodetype=int)
#cs = nx.k_clique_communities(g, 7)
#m = collections.defaultdict(int)
#for c in cs:
#    m[len(c)] += 1
#print m
