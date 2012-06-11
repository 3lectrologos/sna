import numpy
import math
import networkx as nx
import matplotlib.pyplot as plt
import time
import random

class Nodes(nx.MultiGraph):
    P_CONTAM = 0.2
    T_CONTAM = 5
    P_IMMUNITY = 0.3
    P_DEATH = 0.05

    def __init__(self, g):
        nx.MultiGraph.__init__(self, g)
        self.reset()

    def reset(self):
        self.normal = set(self.nodes())
        self.contam = set()
        self.immune = set()
        self.dead = set()
        self.contam_time = {}        

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
        # die
        new_dead = []
        for v in self.contam:
            if random.random() < self.P_DEATH:
                new_dead.append(v)
        [self.contam_time.pop(v) for v in set(new_dead)]
        self.contam = self.contam - set(new_dead)
        self.dead = self.dead | set(new_dead)
        # contaminate
        new_contam = []
        for (v1, v2) in self.edges():
            if v1 in self.contam and v2 in self.normal and \
                    random.random() < self.P_CONTAM:
                new_contam.append(v2)
            elif v2 in self.contam and v1 in self.normal and \
                    random.random() < self.P_CONTAM:
                new_contam.append(v1)
        self.set_contam(new_contam)
        # decontaminate
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
    NORMAL_COLOR = '0.5'
    CONTAM_COLOR = 'g'
    IMMUNE_COLOR = '#E6C200'
    DEAD_COLOR = '0.8'
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

def simulate(g, niter):
#plt.ion()
#draw_network(g, pos)
#time.sleep(1)
    its = []
    dead = []
    immune = []
    normal = []
    it = 0
    for it in range(niter):
        g.reset()
        g.set_contam([1, 42, 50])
        j = 0
        while g.contam:
            g.step()
            j += 1
        its.append(j)
        dead.append(len(g.dead))
        immune.append(len(g.immune))
        normal.append(len(g.normal))
    return {'its': its,
            'dead': dead,
            'immune': immune,
            'normal': normal}
    #draw_network(g, pos)
    #time.sleep(1)

FILE = 'adjnoun.gml'
g = Nodes(nx.read_gml(FILE))
#pos = nx.spring_layout(g, iterations=200)

nsim = 10000
res = simulate(g, nsim)
numnodes = len(g.nodes())

header = FILE + ' (' + str(nsim) + ' simulations)'
print header
print len(header)*'-'
its = res['its']
its_avg = numpy.average(its)
its_std = math.sqrt(numpy.var(its))
dead = res['dead']
dead_avg = 100*numpy.average(dead)/numnodes
dead_std = 100*math.sqrt(numpy.var(dead))/numnodes
immune = res['immune']
immune_avg = 100*numpy.average(immune)/numnodes
immune_std = 100*math.sqrt(numpy.var(immune))/numnodes
normal = res['normal']
normal_avg = 100*numpy.average(normal)/numnodes
normal_std = 100*math.sqrt(numpy.var(normal))/numnodes
print 'avg. its      :', its_avg, '(' + str(its_std) + ')'
print 'avg. dead   % :', dead_avg, '(' + str(dead_std) + ')'
print 'avg. immune % :', immune_avg, '(' + str(immune_std) + ')'
print 'avg. normal % :', normal_avg, '(' + str(normal_std) + ')'

#from pylab import *
#import time
#
#ion()
#
#tstart = time.time()               # for profiling
#x = arange(0,2*pi,0.01)            # x-array
#line, = plot(x,sin(x))
#for i in arange(1,200):
#    line.set_ydata(sin(x+i/10.0))  # update the data
#    draw()                         # redraw the canvas
#
#print 'FPS:' , 200/(time.time()-tstart)
