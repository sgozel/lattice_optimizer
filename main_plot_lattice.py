"""
Copyright 2026 Samuel GOZEL, GNU GPLv3
"""

import os
import re
import math

import lattice_plotter
from lattice_data import BRAVAIS_VECS, SIMU_TORUS

rootdir = os.path.dirname(os.path.abspath(__file__))
latticefilesrootdir = os.path.join(rootdir, 'latticefiles')

def get_filename(lattform, Ns, lattid, typelatt, couplingdetails):
    """
    Extract the full path to a lattice file
    """
    if typelatt == 'original':
        lstr = 'basic_ordering'
        folderpath = os.path.join(latticefilesrootdir, 'original')
    elif typelatt == 'optimized':
        lstr = 'optimized_ordering'
        folderpath = os.path.join(latticefilesrootdir, 'optimized')
    else:
        raise ValueError('Undefined typelatt in get_filename')
    
    pattern = re.compile(rf'{lattform}_{Ns}_{lattid}_{couplingdetails}_{lstr}_bw(\d+)_ops(\d+)\.lattice')
    filefound = False
    for fname in os.listdir(folderpath):
        m = pattern.match(fname)
        if m:
            bw, ops = m.group(1), m.group(2)
            filefound = True
            filename = os.path.join(folderpath, fname)
            break
    if filefound == False:
        print(pattern)
        raise RuntimeError('Lattice file not found')
    return filename, bw, ops

lattices = []

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=16
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'triangle',
    'Ns': int(16),
    'lattid': int(10315),
    'typelatt': 'original', 
    'couplingdetails': 'J',
    'dc': {'J': 1.0}})

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=20
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'triangle',
    'Ns': int(20),
    'lattid': int(10820),
    'typelatt': 'original',
    'couplingdetails': 'J',
    'dc': {'J': 1.0}})

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=21
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'triangle',
    'Ns': int(21),
    'lattid': int(10421),
    'typelatt': 'original',
    'couplingdetails': 'J',
    'dc': {'J': 1.0}})

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=24
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'triangle',
    'Ns': int(24),
    'lattid': int(10424),
    'typelatt': 'optimized',
    'couplingdetails': 'J',
    'dc': {'J': 1.0, 'Jp': math.sqrt(3)}})

#::::::::::::::::::::::::::::::::::::::
# SQUARE NS=25
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'square',
    'Ns': int(25),
    'lattid': int(50005),
    'typelatt': 'original', # 'original'
    'couplingdetails': 'J',
    'dc': {'J': 1.0}})

#::::::::::::::::::::::::::::::::::::::
# TRIANGLE NS=28
#::::::::::::::::::::::::::::::::::::::

lattices.append({
    'lattform': 'triangle',
    'Ns': int(28),
    'lattid': int(20414),
    'typelatt': 'optimized',
    'couplingdetails': 'J',
    'dc': {'J': 1.0, 'Jp': math.sqrt(3)}})

#=============================
# PLOT LATTICES
#=============================

for lattice in lattices:
    
    lattform = lattice['lattform']
    Ns = lattice['Ns']
    lattid = lattice['lattid']
    
    #=============================
    # SEARCH LATTICE FILE
    #=============================
    
    filename, bw, ops = get_filename(lattform,
                                     Ns,
                                     lattid,
                                     lattice['typelatt'],
                                     lattice['couplingdetails'])
    
    #=============================
    # PLOT
    #=============================
    
    latt = lattice_plotter.Lattice(filename, 
                                   dc=lattice['dc'], 
                                   simu_torus=SIMU_TORUS[lattform][lattice['Ns']][lattice['lattid']],
                                   bravais_vecs=BRAVAIS_VECS[lattice['lattform']][lattice['Ns']])
    fig, ax = latt.plot(plot_simu_torus=True,
                        plot_bonds=True,
                        show_numbering=True)
    ax.set_title(f'{lattform.title()}, $N_s = {Ns}$ ({lattid}) [bw={bw}, ops={ops}]')
    # outfigname = os.path.join(rootdir, os.path.basename(filename)[:-8]+'.png')
    # fig.savefig(outfigname, format='png')
