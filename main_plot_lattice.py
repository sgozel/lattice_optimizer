"""
Copyright 2026 Samuel GOZEL, GNU GPLv3
"""

import os
import numpy as np
# import matplotlib.pyplot as plt

import lattice_plotter
from lattice_data import BRAVAIS_VECS, SIMU_TORUS

rootdir = os.path.dirname(os.path.abspath(__file__))
latticefilesrootdir = os.path.join(rootdir, 'latticefiles')

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=16
#::::::::::::::::::::::::::::::::::::::

# lattform = 'triangle'
# Ns = int(16)

# typelatt = 'original'
# file = f'{lattform}_{Ns}_J_basic_ordering_bw13_ops384'
# dc = {'J': 1.0}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_J_optimized_ordering_bw12_ops360'
# dc = {'J': 1.0}

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=20
#::::::::::::::::::::::::::::::::::::::

# lattform = 'triangle'
# Ns = int(20)

# typelatt = 'original'
# file = f'{lattform}_{Ns}_J_basic_ordering_bw19_ops484'
# dc = {'J': 1.0}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_J_optimized_ordering_bw9_ops460'
# dc = {'J': 1.0}

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=24
#::::::::::::::::::::::::::::::::::::::

lattform = 'triangle'
Ns = int(24)

# typelatt = 'original'
# file = f'{lattform}_{Ns}_J_basic_ordering_bw23_ops904'
# dc = {'J': 1.0}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_J_optimized_ordering_bw13_ops672'
# dc = {'J': 1.0}

# typelatt = 'original'
# file = f'{lattform}_24_JJp_basic_ordering_bw23_ops2088'
# dc = {'J': 1.0, 'Jp': np.sqrt(3)}

typelatt = 'optimized'
file = f'{lattform}_24_JJp_optimized_ordering_bw17_ops1784'
dc = {'J': 1.0, 'Jp': np.sqrt(3)}

#::::::::::::::::::::::::::::::::::::::
# SQUARE NS=25
#::::::::::::::::::::::::::::::::::::::

# lattform = 'square'
# Ns = int(25)

# typelatt = 'original'
# file = f'{lattform}_{Ns}_J_basic_ordering_bw20_ops430'
# dc = {'J': 1.0}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_J_optimized_ordering_bw14_ops414'
# dc = {'J': 1.0}

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=28
#::::::::::::::::::::::::::::::::::::::

# lattform = 'triangle'
# Ns = int(28)

# typelatt = 'original'
# file = f'{lattform}_{Ns}_J_basic_ordering_bw27_ops1076'
# dc = {'J': 1.0}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_J_optimized_ordering_bw18_ops884'
# dc = {'J': 1.0}

# typelatt = 'original'
# file = f'{lattform}_{Ns}_JJp_basic_ordering_bw27_ops2596'
# dc = {'J': 1.0, 'Jp': np.sqrt(3)}

# typelatt = 'optimized'
# file = f'{lattform}_{Ns}_JJp_optimized_ordering_bw22_ops2368'
# dc = {'J': 1.0, 'Jp': np.sqrt(3)}

#=============================
# PLOT
#=============================

filename_fullpath = os.path.join(latticefilesrootdir, typelatt, 
                                 file + '.lattice')
latt = lattice_plotter.Lattice(filename_fullpath, 
                               dc=dc, 
                               simu_torus=SIMU_TORUS[lattform][Ns],
                               bravais_vecs=BRAVAIS_VECS[lattform][Ns])
fig, ax = latt.plot(plot_simu_torus=True,
                    plot_bonds=True,
                    show_numbering=True)
ax.set_title(f'{lattform.title()}, $N_s = {Ns}$')
# outfigname = os.path.join(rootdir, file+'.png')
# fig.savefig(outfigname, format='png')
