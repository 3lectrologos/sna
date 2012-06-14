"""
ETH Zurich | Spring Semester 2012 | Introduction to Social Network Analysis
Final Project: Simulating disease spread in social networks
"""

__author__ = "Alkis Gkotovos"
__email__ = "alkisg@student.ethz.ch"


import numpy
import math
import networkx as nx
import matplotlib.pyplot as plt
import time
import random
import os
import re
import operator
import progressbar as pb

NORMAL_COLOR = '#5555EE'
CONTAM_COLOR = '#009926'
IMMUNE_COLOR = '#E6C200'
DEAD_COLOR = '#BBBBBB'
FONT_SIZE = 12
FIG_PATH = os.getcwd() + '/report/figures/'

def clamp(val, minimum=0, maximum=255):
    if val < minimum:
        return minimum
    if val > maximum:
        return maximum
    return val

def colorscale(hexstr, scalefactor):
    """
    Scale a hex string by ``scalefactor``.
    """
    hexstr = hexstr.strip('#')
    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr
    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)
    r = clamp(r * scalefactor)
    g = clamp(g * scalefactor)
    b = clamp(b * scalefactor)
    return "#%02x%02x%02x" % (r, g, b)

class Simr(nx.Graph):
    """
    Encapsulate a graph to allow for disease spread simulation.
    """
    # Probability of contaminating a healthy neighbor
    P_CONTAM = 0.05
    # Time-steps for which the disease is active
    T_CONTAM = 6
    # Probability that a node will become immune after the disease has passed
    P_IMMUNITY = 0.2
    # Probability that a node will die at each time step that it is diseased
    P_DEATH = 0.05
    # Change this for starting with given immune nodes
    init_immune = []

    # Initialize with the given graph
    def __init__(self, g, k=0.05, pos=None):
        nx.Graph.__init__(self, g)
        self.gtype = 'UNK'
        self.gparam = []
        self.pos = pos
        self.k = k
        self.reset()

    # Initialize with a generated graph of three classes
    def __init__(self, gtype='ER', gparam=[100, 0.1], k=0.05, pos=None):
        if gtype == 'ER':
            g = nx.erdos_renyi_graph(*gparam)
        elif gtype == 'WS':
            g = nx.watts_strogatz_graph(*gparam)
        elif gtype == 'BA':
            g = nx.barabasi_albert_graph(*gparam)
        else:
            raise Exception("Unknown graph model type.")
        nx.Graph.__init__(self, g)
        self.gtype = gtype
        self.gparam = gparam
        self.pos = pos
        self.k = k
        self.reset()

    def reset(self):
        self.normal = set(self.nodes())
        self.contam = set()
        self.immune = set()
        self.dead = set()
        self.contam_time = {}
        n = len(self.nodes())
        self.set_contam(random.sample(self.nodes(),
                                      int(math.ceil(self.k*n))))
        self.set_immune(self.init_immune)

    def set_normal(self, new):
        new = set(new)
        self.normal = self.normal | new
        self.contam = self.contam - new
        self.immune = self.immune - new
        self.dead = self.dead - new

    def set_contam(self, new):
        new = set(new)
        self.normal = self.normal - new
        self.contam = self.contam | new
        for v in new:
            self.contam_time[v] = self.T_CONTAM
        self.immune = self.immune - new
        self.dead = self.dead - new

    def set_immune(self, new):
        new = set(new)
        self.normal = self.normal - new
        self.contam = self.contam - new
        self.immune = self.immune | new
        self.dead = self.dead - new

    def set_dead(self, new):
        new = set(new)
        self.normal = self.normal - new
        self.contam = self.contam - new
        self.immune = self.immune - new
        self.dead = self.dead | new

    def nondead_edges(self):
        return [v for v in self.edges() if
                v[0] not in self.dead and
                v[1] not in self.dead]

    def step(self):
        # For every contaminated node randomly decide if he dies
        # according to P_DEATH.
        new_dead = []
        for v in self.contam:
            if random.random() < self.P_DEATH:
                new_dead.append(v)
        [self.contam_time.pop(v) for v in set(new_dead)]
        self.contam = self.contam - set(new_dead)
        self.dead = self.dead | set(new_dead)
        # For every neighbor of a contaminated node randomly decide if he
        # catches the disease according to P_CONTAM.
        new_contam = []
        for (v1, v2) in self.edges():
            if v1 in self.contam and v2 in self.normal and \
                    random.random() < self.P_CONTAM:
                new_contam.append(v2)
            elif v2 in self.contam and v1 in self.normal and \
                    random.random() < self.P_CONTAM:
                new_contam.append(v1)
        self.set_contam(new_contam)
        # For every contaminated node whose contamination period is over
        # (according to T_CONTAM), randomly decide if he becomes immune
        # or just healthy according to P_IMMUNITY.
        new_immune = []
        for (v, t) in self.contam_time.iteritems():
            if t == 0:
                new_immune.append(v)
            else:
                self.contam_time[v] -= 1
        [self.contam_time.pop(v) for v in set(new_immune)]
        self.contam = self.contam - set(new_immune)
        for v in new_immune:
            if random.random() < self.P_IMMUNITY:
                self.set_immune([v])
            else:
                self.set_normal([v])

    def plot(self):
        if self.pos == None:
            self.pos = nx.graphviz_layout(self)
        NODE_SIZE = 500
        plt.clf()
        nx.draw_networkx_nodes(self, pos=self.pos,
                               nodelist=self.normal,
                               node_color=NORMAL_COLOR,
                               node_size=NODE_SIZE)
        nx.draw_networkx_nodes(self, pos=self.pos,
                               nodelist=self.contam,
                               node_color=CONTAM_COLOR,
                               node_size=NODE_SIZE)
        nx.draw_networkx_nodes(self, pos=self.pos,
                               nodelist=self.immune,
                               node_color=IMMUNE_COLOR,
                               node_size=NODE_SIZE)
        nx.draw_networkx_nodes(self, pos=self.pos,
                               nodelist=self.dead,
                               node_color=DEAD_COLOR,
                               node_size=NODE_SIZE)
        nx.draw_networkx_edges(self, pos=self.pos,
                               edgelist=self.nondead_edges(),
                               width=2,
                               edge_color='0.2')
        nx.draw_networkx_labels(self, pos=self.pos,
                                font_color='0.95', font_size=11)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.draw()

    def print_summary(self):
        print 'normal :', len(self.normal)
        print 'contam :', len(self.contam)
        print 'immune :', len(self.immune)
        print 'dead   :', len(self.dead)

    def graph_type(self):
        res = self.gtype + '('
        for i in range(len(self.gparam) - 1):
            res += str(self.gparam[i]) + ', '
        res += str(self.gparam[-1]) + ')'
        return res

    def graph_name(self):
        res = self.gtype + '_'
        for i in range(len(self.gparam) - 1):
            res += str(self.gparam[i]) + '_'
        res += str(self.gparam[-1])
        res = re.sub('\.', '', res)
        return res

    def print_graph_info(self):
        ac = nx.average_clustering(self)
        ap = nx.average_shortest_path_length(self)
        print self.graph_type()
        print 'Avg. clustering =', ac
        print 'Avg. short. path len. =', ap

    def __str__(self):
        ret = ''
        ret += 'normal : ' + str(list(self.normal)) + '\n'
        ret += 'contam : ' + str(list(self.contam)) + '\n'
        ret += 'immune : ' + str(list(self.immune)) + '\n'
        ret += 'dead   : ' + str(list(self.dead)) + '\n'
        return ret

def simulate_evo(g, save=False, savename='graph_evo', show=False, freq=10**8):
    """
    Simulate a single evolution of the disease in the given network.
    The evolution ends when there are no more diseased nodes. Return
    information about the network at each time step.
    """
    its = 0
    contam = [len(g.contam)]
    dead = [len(g.dead)]
    immune = [len(g.immune)]
    normal = [len(g.normal)]
    g.reset()
    if save or show:
        g.plot()
        if save:
            savefig(savename + '_init')
        if show:
            plt.show()
    while g.contam:
        g.step()
        contam.append(len(g.contam))
        dead.append(len(g.dead))
        immune.append(len(g.immune))
        normal.append(len(g.normal))
        its += 1
        if its % freq == 0 and (save or show):
            g.plot()
            if save:
                savefig(savename + '_' + str(its))
            if show:
                plt.show()
    if save or show:
        g.plot()
        if save:
            savefig(savename + '_final')
        if show:
            plt.show()
    return {'its': its,
            'nodes': len(g.nodes()),
            'contam': contam,
            'dead': dead,
            'immune': immune,
            'normal': normal}

def simulate(g, niter, plot=False):
    """
    Simulate niter evolutions of the disease in the given network.
    Each evolution ends when there are no more diseased nodes. Return
    information about the final state of the network at each iteration.
    """
    its = []
    dead = []
    immune = []
    max_contam = []
    normal = []
    for it in range(niter):
        g.reset()
        mc = 0
        j = 0
        while g.contam:
            g.step()
            mc = max(mc, len(g.contam))
            j += 1
            if plot:
                g.plot()
                plt.show()
        its.append(j)
        dead.append(len(g.dead))
        immune.append(len(g.immune))
        max_contam.append(mc)
        normal.append(len(g.normal))
    return {'its': its,
            'max_contam': max_contam,
            'nodes': len(g.nodes()),
            'dead': dead,
            'immune': immune,
            'normal': normal}

def plot_evo(data):
    """
    Create line plot using data returned by ``simulate_evo``.
    """
    its = data['its']
    nodes = data['nodes']
    contam = [(100.0*d)/nodes for d in data['contam']]
    dead = [(100.0*d)/nodes for d in data['dead']]
    immune = [(100.0*d)/nodes for d in data['immune']]
    normal = [(100.0*d)/nodes for d in data['normal']]
    legitems = []
    plt.clf()
    p = plt.plot(normal, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(NORMAL_COLOR, 0.4),
                 markersize=6,
                 color=NORMAL_COLOR)
    legitems.append(p)
    p = plt.plot(contam, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(CONTAM_COLOR, 0.2),
                 markersize=6,
                 color=CONTAM_COLOR)
    legitems.append(p)
    p = plt.plot(dead, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(DEAD_COLOR, 0.1),
                 markersize=6,
                 color=colorscale(DEAD_COLOR, 0.7))
    legitems.append(p)
    p = plt.plot(immune, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(IMMUNE_COLOR, 0.4),
                 markersize=6,
                 color=IMMUNE_COLOR)
    legitems.append(p)
    plt.xlim([0, its])
    plt.ylim([0, 100])
    plt.xlabel('Time step')
    plt.ylabel('Node percentage')
    add_legend(legitems,
               ['healthy', 'infected', 'removed', 'immune'],
               1)

def average_data(data):
    """
    Find mean and std. deviation of data returned by ``simulate``.
    """
    numnodes = data['nodes']
    its = data['its']
    its_mean = numpy.average(its)
    its_std = math.sqrt(numpy.var(its))
    dead = data['dead']
    dead_mean = 100.0*numpy.average(dead)/numnodes
    dead_std = 100.0*math.sqrt(numpy.var(dead))/numnodes
    immune = data['immune']
    immune_mean = 100.0*numpy.average(immune)/numnodes
    immune_std = 100.0*math.sqrt(numpy.var(immune))/numnodes
    max_contam = data['max_contam']
    max_contam_mean = 100.0*numpy.average(max_contam)/numnodes
    max_contam_std = 100.0*math.sqrt(numpy.var(max_contam))/numnodes
    normal = data['normal']
    normal_mean = 100.0*numpy.average(normal)/numnodes
    normal_std = 100.0*math.sqrt(numpy.var(normal))/numnodes
    return {'its': (its_mean, its_std),
            'nodes': numnodes,
            'dead': (dead_mean, dead_std),
            'immune': (immune_mean, immune_std),
            'max_contam': (max_contam_mean, max_contam_std),
            'normal': (normal_mean, normal_std)}

def simulate_range(gfun, vals, niter):
    """
    Simulate disease evolution for a range of different values of a parameter.
    """
    data = []
    widgets = [pb.Bar('>'), ' ', pb.ETA(), ' ', pb.ReverseBar('<')]
    progress = pb.ProgressBar(widgets=widgets)
    for val in progress(vals):
        g = gfun(val)
        r = simulate(g, niter)
        data.append(average_data(r))
    return {'vals': vals, 'data': data}

def plot_range(data, xlabel='x', ylabel='Average node percentage'):
    """
    Plot network evolution data returned by ``simulate_range``.
    """
    vals = data['vals']
    data = data['data']
    dead = [r['dead'][0] for r in data]
    dead_std = [r['dead'][1] for r in data]
    immune = [r['immune'][0] for r in data]
    immune_std = [r['immune'][1] for r in data]
    max_contam = [r['max_contam'][0] for r in data]
    max_contam_std = [r['max_contam'][1] for r in data]
    legitems = []
    plt.clf()
    p = plt.errorbar(vals, max_contam, yerr=max_contam_std,
                 linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(CONTAM_COLOR, 0.2),
                 markersize=6,
                 color=CONTAM_COLOR)
    legitems.append(p)
    p = plt.errorbar(vals, dead, yerr=dead_std,
                     linestyle='-', marker='o', linewidth=2,
                     markeredgecolor=colorscale(DEAD_COLOR, 0.1),
                     markersize=6,
                     color=colorscale(DEAD_COLOR, 0.7))
    legitems.append(p)
    p = plt.errorbar(vals, immune, yerr=immune_std,
                 linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(IMMUNE_COLOR, 0.4),
                 markersize=6,
                 color=IMMUNE_COLOR)
    legitems.append(p)
    plt.xlim([0, vals[-1]])
    plt.ylim([0, 100])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    add_legend(legitems, ['max. infected', 'dead', 'immune'], 4)

def add_legend(legitems, strings, pos):
    plt.legend(legitems, strings, pos)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    llines = leg.get_lines()
    frame  = leg.get_frame()
    frame.set_facecolor('1.0')
    plt.setp(ltext, fontsize=FONT_SIZE)
    plt.setp(llines, linewidth=1.0)

def savefig(outfile):
    """
    Save current figure to file.
    """
    # General font size
    font = {'family': 'sans', 'weight': 'normal', 'size': FONT_SIZE}
    plt.rc('font', **font)
    # Adjust size
    plt.gcf().set_size_inches(8, 6)
    # Save
    outfile = FIG_PATH + outfile
    plt.savefig(outfile + '.pdf', format='pdf', bbox_inches='tight') 

def gen_er_plots():
    # ER graph state evolution plots
    g = Simr(gtype='ER', gparam=[50, 0.3], k=0.05)
    g.k = 0.05
    res = simulate_evo(g, save=True, savename='evo_' + g.graph_name(), freq=10)
    plt.close()

    # ER evolution plots
    g3 = Simr(gtype='ER', gparam=[500, 0.015], k=0.05)
    g3.print_graph_info()
    res = simulate_evo(g3)
    plot_evo(res)
    savefig('evo_' + g3.graph_name())
    g4 = Simr(gtype='ER', gparam=[500, 0.03], k=0.05)
    g4.print_graph_info()
    res = simulate_evo(g4)
    plot_evo(res)
    savefig('evo_' + g4.graph_name())

    # ER summary plots
    NITER = 2
    gfun = lambda p: Simr(gtype='ER', gparam=[50, p], k=0.05)
    res = simulate_range(gfun, numpy.linspace(0, 0.8, 20), NITER)
    plot_range(res, xlabel='p')
    savefig('sum_ER_50_p')

    gfun = lambda p: Simr(gtype='ER', gparam=[500, p], k=0.05)
    res = simulate_range(gfun, numpy.linspace(0, 0.2, 20), NITER)
    plot_range(res, xlabel='p')
    savefig('sum_ER_500_p')

def gen_ws_plots():
    # WS evolution plots
    g3 = Simr(gtype='WS', gparam=[500, 15, 0.1], k=0.05)
    g3.print_graph_info()
    res = simulate_evo(g3)
    plot_evo(res)
    savefig('evo_' + g3.graph_name())
    g4 = Simr(gtype='WS', gparam=[500, 30, 0.1], k=0.05)
    g4.print_graph_info()
    res = simulate_evo(g4)
    plot_evo(res)
    savefig('evo_' + g4.graph_name())

    # WS summary plots
    NITER = 2
    gfun = lambda p: Simr(gtype='WS', gparam=[500, 30, p], k=0.05)
    res = simulate_range(gfun, numpy.linspace(0, 1, 20), NITER)
    plot_range(res, xlabel='p')
    savefig('sum_WS_500_25_p')
    gfun = lambda k: Simr(gtype='WS', gparam=[500, k, 0.1], k=0.05)
    res = simulate_range(gfun, range(3, 40), NITER)
    plot_range(res, xlabel='k')
    savefig('sum_WS_500_k_01')

def gen_ba_plots():
    # BA evolution plots
    g = Simr(gtype='BA', gparam=[500, 3], k=0.05)
    g.print_graph_info()
    res = simulate_evo(g)
    plot_evo(res)
    savefig('evo_' + g.graph_name())
    g = Simr(gtype='BA', gparam=[500, 7], k=0.05)
    g.print_graph_info()
    res = simulate_evo(g)
    plot_evo(res)
    savefig('evo_' + g.graph_name())

    # BA summary plots
    NITER = 2
    gfun = lambda p: Simr(gtype='BA', gparam=[500, p], k=0.05)
    res = simulate_range(gfun, range(1, 25), NITER)
    plot_range(res, xlabel='m')
    savefig('sum_BA_500_m')

def rm_hubs(g, n):
    """
    Make ``n'' nodes with the largest closeness centrality immune.
    """
    c = nx.closeness_centrality(g)
    s = sorted(c.iteritems(), key=operator.itemgetter(1), reverse=True)
    hubs = [h[0] for h in s[:n]]
    g.init_immune = hubs
    return g

def gen_rmhubs_plots():
    NITER = 10
    # ER
    gfun = lambda n: rm_hubs(Simr(gtype='ER', gparam=[500, 0.05], k=0.05), n)
    res = simulate_range(gfun, range(0, 31, 3), NITER)
    plot_range(res, xlabel='Number of removed nodes')
    savefig('hubs_ER_500_005')

    # WS
    gfun = lambda n: rm_hubs(Simr(gtype='WS', gparam=[500, 29, 0.1], k=0.05), n)
    res = simulate_range(gfun, range(0, 31, 3), NITER)
    plot_range(res, xlabel='Number of removed nodes')
    savefig('hubs_WS_500_29_01')

    # BA
    gfun = lambda n: rm_hubs(Simr(gtype='BA', gparam=[500, 13], k=0.05), n)
    res = simulate_range(gfun, range(0, 31, 3), NITER)
    plot_range(res, xlabel='Number of removed nodes')
    savefig('hubs_BA_500_13')

if __name__=="__main__":
    # ER plots
    gen_er_plots()
    # WS plots
    gen_ws_plots()
    # BA plots
    gen_ba_plots()
    # Removed hubs plots
    gen_rmhubs_plots()
