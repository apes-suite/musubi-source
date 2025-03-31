import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 20})


# Read data
D_p = 15.0e-3 
clrs = ['#427299', '#009aa3','#4bb75f','#e7b70a']

class ParticleData:
    def __init__(self, infile, label):
        self.infile = infile
        self.data = self.read_data(self.infile)
        self.label = label

    def read_data(self,infile):
        d = pd.read_csv(infile, skiprows = 1, 
                    sep='\s+')
        return d

    def plot_trajectory(self, fig, ax, line_style='-', clr='k'):
        ax.plot(self.data['t'], self.data['z']/D_p, linestyle=line_style, 
                color=clr, label=self.label)


class ExperimentalParticleData(ParticleData):
    def read_data(self,infile):
        d = pd.read_csv(infile, header=None, 
                    names=['t','z'],
                    sep=', ')
        # Add 0.5 because we're comparing the CENTER of the sphere 
        d['z'] = (d['z'] + 0.5)*D_p   
        return d

    def plot_trajectory(self, fig, ax, line_style='-', clr='k'):
        ax.plot(self.data['t'], self.data['z']/D_p, 
                linestyle=line_style, color=clr, marker='s',
                label=self.label)

# Data from fully-resolved simulations using Discrete Particle Simulations (DPS)
DPS_Re11_6 = ParticleData('./trajectories/particle0000001.dat', 'DPS')


# Experimental data 
# Taken from the plots in 
#
# [1] A. ten Cate, C. H. Nieuwstad, J. J. Derksen, and H. E. A. Van den Akker, 
# “Particle imaging velocimetry experiments and lattice-Boltzmann simulations 
# on a single sphere settling under gravity,” 
# Physics of Fluids, vol. 14, no. 11, pp. 4012–4025, Nov. 2002, doi: 10.1063/1.1512918.
#
# using WebPlotDigitizer 
EXP_11_6 = ExperimentalParticleData('./reference_data/tencate_Re11_6.csv', 'Ten Cate')

fig, ax = plt.subplots()
DPS_Re11_6.plot_trajectory(fig, ax, '-',clrs[0])
EXP_11_6.plot_trajectory(fig, ax,'',clrs[0])

ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$z/D_p$')
ax.grid()
ax.legend()
fig.set_size_inches((10,8))
plt.tight_layout()
plt.savefig(f'tencate_DPS_RE12_EXP.png')

