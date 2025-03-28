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
                    delim_whitespace=True)
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

# Data from fully-resolved simulations using Momentum-Exchange-Method (MEM)
MEM_Re32_2 = ParticleData('../Re32_2/trajectories/particle0001.dat', '_MEM Re=32.2')
MEM_Re11_6 = ParticleData('../Re11_6/trajectories/particle0001.dat', '_MEM Re=11.6')
MEM_Re4_1 = ParticleData('../Re4_1/trajectories/particle0001.dat', '_MEM Re=4.1')
MEM_Re1_5 = ParticleData('../Re1_5/trajectories/particle0001.dat', '_MEM Re=1.5')

# Data from unresolved simulations using Discrete Partice Simulations (DPS)
DPS_32_2 = ParticleData('/Users/vlogmantg/Documents/PhD/Code/mus_simulations/tencate_DPS/Re32_2/trajectories/particle0000001.dat', '_DPS Re=32.2')
DPS_11_6 = ParticleData('/Users/vlogmantg/Documents/PhD/Code/mus_simulations/tencate_DPS/Re11_6/trajectories/particle0000001.dat', '_DPS Re=11.6')
DPS_4_1 = ParticleData('/Users/vlogmantg/Documents/PhD/Code/mus_simulations/tencate_DPS/Re4_1/trajectories/particle0000001.dat', '_DPS Re=4.1')
DPS_1_5 = ParticleData('/Users/vlogmantg/Documents/PhD/Code/mus_simulations/tencate_DPS/Re1_5/trajectories/particle0000001.dat', '_DPS Re=1.5')

# Experimental data 
# Taken from the plots in 
#
# [1] A. ten Cate, C. H. Nieuwstad, J. J. Derksen, and H. E. A. Van den Akker, 
# “Particle imaging velocimetry experiments and lattice-Boltzmann simulations 
# on a single sphere settling under gravity,” 
# Physics of Fluids, vol. 14, no. 11, pp. 4012–4025, Nov. 2002, doi: 10.1063/1.1512918.
#
# using WebPlotDigitizer 
EXP_32_2 = ExperimentalParticleData('../experimental_data/tencate_Re32_2.csv', 'Re=32.2')
EXP_11_6 = ExperimentalParticleData('../experimental_data/tencate_Re11_6.csv', 'Re=11.6')
EXP_4_1 = ExperimentalParticleData('../experimental_data/tencate_Re4_1.csv', 'Re=4.1')
EXP_1_5 = ExperimentalParticleData('../experimental_data/tencate_Re1_5.csv', 'Re=1.5')

fig, ax = plt.subplots()
MEM_Re32_2.plot_trajectory(fig, ax, '-',clrs[0])
MEM_Re11_6.plot_trajectory(fig, ax, '-',clrs[1])
MEM_Re4_1.plot_trajectory(fig, ax, '-',clrs[2])
MEM_Re1_5.plot_trajectory(fig, ax, '-',clrs[3])

DPS_32_2.plot_trajectory(fig, ax, '--',clrs[0])
DPS_11_6.plot_trajectory(fig, ax, '--',clrs[1])
DPS_4_1.plot_trajectory(fig, ax,'--',clrs[2])
DPS_1_5.plot_trajectory(fig, ax,'--',clrs[3])

EXP_32_2.plot_trajectory(fig, ax,'',clrs[0])
EXP_11_6.plot_trajectory(fig, ax,'',clrs[1])
EXP_4_1.plot_trajectory(fig, ax,'',clrs[2])
EXP_1_5.plot_trajectory(fig, ax,'',clrs[3])

ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$z/D_p$')
ax.grid()
ax.legend()
fig.set_size_inches((10,8))
plt.tight_layout()
plt.savefig('tencate_MEM_DPS_EXP.pdf')
plt.show()

