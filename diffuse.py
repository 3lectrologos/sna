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
import progressbar as pb

NORMAL_COLOR = '#5555EE'
CONTAM_COLOR = '#009926'
IMMUNE_COLOR = '#E6C200'
DEAD_COLOR = '#DDDDDD'

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

class Nodes(nx.MultiGraph):
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

    def __init__(self, g):
        nx.MultiGraph.__init__(self, g)
        self.reset()

    def __init__(self, g, k):
        nx.MultiGraph.__init__(self, g)
        self.k = k
        self.reset()

    def reset(self):
        self.normal = set(self.nodes())
        self.contam = set()
        self.immune = set()
        self.dead = set()
        self.contam_time = {}
        self.set_contam(random.sample(self.nodes(), self.k))

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

    def print_summary(self):
        print 'normal :', len(self.normal)
        print 'contam :', len(self.contam)
        print 'immune :', len(self.immune)
        print 'dead   :', len(self.dead)

    def __str__(self):
        ret = ''
        ret += 'normal : ' + str(list(self.normal)) + '\n'
        ret += 'contam : ' + str(list(self.contam)) + '\n'
        ret += 'immune : ' + str(list(self.immune)) + '\n'
        ret += 'dead   : ' + str(list(self.dead)) + '\n'
        return ret

def draw_network(g, pos):
    NODE_SIZE = 500
    plt.clf()
    nx.draw_networkx_nodes(g, pos=pos,
                           nodelist=g.normal,
                           node_color=NORMAL_COLOR,
                           node_size=NODE_SIZE)
    nx.draw_networkx_nodes(g, pos=pos,
                           nodelist=g.contam,
                           node_color=CONTAM_COLOR,
                           node_size=NODE_SIZE)
    nx.draw_networkx_nodes(g, pos=pos,
                           nodelist=g.immune,
                           node_color=IMMUNE_COLOR,
                           node_size=NODE_SIZE)
    nx.draw_networkx_nodes(g, pos=pos,
                           nodelist=g.dead,
                           node_color=DEAD_COLOR,
                           node_size=NODE_SIZE)
    nx.draw_networkx_edges(g, pos=pos,
                           edgelist=g.nondead_edges(),
                           width=2,
                           edge_color='0.2')
    nx.draw_networkx_labels(g, pos=pos, font_color='0.95', font_size=11)
    plt.draw()

def simulate_evo(g, plot=False):
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
    while g.contam:
        g.step()
        contam.append(len(g.contam))
        dead.append(len(g.dead))
        immune.append(len(g.immune))
        normal.append(len(g.normal))
        its += 1
        if plot:
            draw_network(g, pos)
            time.sleep(1)
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
    it = 0
    widgets = [pb.Bar('>'), ' ', pb.ETA(), ' ', pb.ReverseBar('<')]
    progress = pb.ProgressBar(widgets=widgets)
    for it in progress(range(niter)):
        g.reset()
        mc = 0
        j = 0
        while g.contam:
            g.step()
            mc = max(mc, len(g.contam))
            j += 1
            if plot:
                draw_network(g, pos)
                time.sleep(1)
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
    Create line plot depicting a single network evolution.
    """
    its = data['its']
    nodes = data['nodes']
    contam = [(1.0*d)/nodes for d in data['contam']]
    dead = [(1.0*d)/nodes for d in data['dead']]
    immune = [(1.0*d)/nodes for d in data['immune']]
    normal = [(1.0*d)/nodes for d in data['normal']]
    legitems = []
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
    p = plt.plot(normal, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(NORMAL_COLOR, 0.4),
                 markersize=6,
                 color=NORMAL_COLOR)
    legitems.append(p)
    plt.xlim([0, its])
    plt.ylim([0, 1])
    plt.legend(legitems, ['contaminated %', 'dead %', 'immune %', 'normal %'], 1)
    plt.show()

def average_data(data):
    """
    Find mean and std. deviation of data returned by ``simulate``.
    """
    numnodes = data['nodes']
    its = data['its']
    its_mean = numpy.average(its)
    its_std = math.sqrt(numpy.var(its))
    dead = data['dead']
    dead_mean = 1.0*numpy.average(dead)/numnodes
    dead_std = 1.0*math.sqrt(numpy.var(dead))/numnodes
    immune = data['immune']
    immune_mean = 1.0*numpy.average(immune)/numnodes
    immune_std = 1.0*math.sqrt(numpy.var(immune))/numnodes
    max_contam = data['max_contam']
    max_contam_mean = 1.0*numpy.average(max_contam)/numnodes
    max_contam_std = 1.0*math.sqrt(numpy.var(max_contam))/numnodes
    normal = data['normal']
    normal_mean = 1.0*numpy.average(normal)/numnodes
    normal_std = 1.0*math.sqrt(numpy.var(normal))/numnodes
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
    for val in vals:
        g = gfun(val)
        r = simulate(g, niter)
        data.append(average_data(r))
    return {'vals': vals, 'data': data}

def plot_range(data, xlabel='x', ylabel='y'):
    """
    Plot network evolution data returned by ``simulate_range``.
    """
    vals = data['vals']
    data = data['data']
    dead = [r['dead'][0] for r in data]
    immune = [r['immune'][0] for r in data]
    max_contam = [r['max_contam'][0] for r in data]
    legitems = []
    p = plt.plot(vals, max_contam, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(CONTAM_COLOR, 0.2),
                 markersize=6,
                 color=CONTAM_COLOR)
    legitems.append(p)
    p = plt.plot(vals, dead, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(DEAD_COLOR, 0.1),
                 markersize=6,
                 color=colorscale(DEAD_COLOR, 0.7))
    legitems.append(p)
    print dead
    p = plt.plot(vals, immune, linestyle='-', marker='o', linewidth=2,
                 markeredgecolor=colorscale(IMMUNE_COLOR, 0.4),
                 markersize=6,
                 color=IMMUNE_COLOR)
    legitems.append(p)
    #p = plt.plot(normal, linestyle='-', marker='o', linewidth=2,
    #             markeredgecolor=colorscale(NORMAL_COLOR, 0.4),
    #             markersize=6,
    #             color=NORMAL_COLOR)
    #legitems.append(p)
    #plt.xlim([0, vals])
    plt.ylim([0, 1])
    plt.legend(legitems, ['max. contam. %', 'dead %', 'immune %'], 1)
    plt.show()

if __name__=="__main__":
    gfun = lambda n: Nodes(nx.barabasi_albert_graph(100, n), k=3)
    res = simulate_range(gfun, [2, 3, 4, 5, 10, 15, 20, 25, 30], 10000)
    plot_range(res)
    #g = Nodes(nx.erdos_renyi_graph(100, 0.9))
    #g = Nodes(nx.watts_strogatz_graph(100, 30, 0.8))
    #g = Nodes(nx.barabasi_albert_graph(10000, 5))
    #pos = nx.spring_layout(g, iterations=100)
    #plt.ion()

    #res = simulate_evo(g, plot=False)
    #plot_evo(res)

    #nsim = 1000
    #res = simulate(g, nsim, plot=True)
    #numnodes = len(g.nodes())
    
    #NAME = 'Erdos-Renyi (100, 0.3)'
    #header = NAME + ' (' + str(nsim) + ' simulations)'
    #print header
    #print len(header)*'-'
 
    #print 'avg. its      :', its_avg, '(' + str(its_std) + ')'
    #print 'avg. dead   % :', dead_avg, '(' + str(dead_std) + ')'
    #print 'avg. immune % :', immune_avg, '(' + str(immune_std) + ')'
    #print 'avg. normal % :', normal_avg, '(' + str(normal_std) + ')'
