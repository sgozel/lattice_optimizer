"""
Copyright 2026 Samuel GOZEL, GNU GPLv3
"""

import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt

   
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern",
    "font.size": 16,
    'axes.labelsize': 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 12,
})


class Lattice:
    """
    Class to store and plot a lattice - a set of sites and a set of bonds
    
    """
    
    def __init__(self, filepath, dc, simu_torus, bravais_vecs):
        """
        Constructor

        Parameters
        ----------
        filepath : str
            Path to lattice file
        dc : dict
            Bond distance by coupling
        simu_torus : dict
            Simulation torus - dictonnary with 2 keys: 't1' and 't2'values, with
            values being numpy arrays expressing the simulation torus in the {a1, a2}
            Bravais basis
        bravais_vecs : dict
            Bravais basis vectors - dictionnary with 2 keys: 'a1' and 'a2', with
            values being numpy arrays expressing the Bravais basis vectors in the
            canonical basis
        
        """
        if not os.path.exists(filepath):
            raise ValueError(f'Lattice file {filepath} does not exist.')
        
        self.filepath = filepath
        self.simu_torus = simu_torus
        self.bravais_vecs = bravais_vecs
        self.dc = dc
        
        # Simulation torus in the canonical basis
        self.t1 = self.simu_torus['t1'][0] * self.bravais_vecs['a1'] + self.simu_torus['t1'][1] * self.bravais_vecs['a2']
        self.t2 = self.simu_torus['t2'][0] * self.bravais_vecs['a1'] + self.simu_torus['t2'][1] * self.bravais_vecs['a2']
        
        # read lattice file
        self._read_lines()
        self._parse()
        
        self.sites_x = []
        self.sites_y = []
        for i in self.sites.keys():
            self.sites_x.append(self.sites[i][0])
            self.sites_y.append(self.sites[i][1])
        self.sites_x = np.array(self.sites_x)
        self.sites_y = np.array(self.sites_y)
        self.minX = np.min(self.sites_x)
        self.maxX = np.max(self.sites_x)
        self.minY = np.min(self.sites_y)
        self.maxY = np.max(self.sites_y)
        
        self.ghosts = []
        self._get_ghosts()
        
        return
    
    
    def _read_lines(self):
        """
        Read the lattice file

        """
        with open(self.filepath, encoding="utf-8", errors="replace") as f:
            self.lines = [line.rstrip("\n") for line in f]
        return
    
    
    def _parse(self):
        """
        Parse lines from the lattice file
        
        """
        #============================
        # Parse number of sites
        #============================
        pattern = r"\[sites\]=(\d+)"
        idx = int(0)
        m = re.match(pattern, self.lines[idx])
        if not m:
            raise ValueError(f'Expected [sites]=<n> on first line. Got: {self.lines[idx]!r}')
        self.Ns = int(m.group(1))
        
        #============================
        # Parse sites coordinates
        #============================
        idx += 1
        self.sites = {}
        for _ in range(self.Ns):
            parts = self.lines[idx].split()
            if len(parts) != 4:
                raise ValueError(f"Expected '<index> <x> <y> <z>', got: {self.lines[idx]!r}")
            site_idx = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            self.sites[site_idx] = np.array([x, y, z])
            idx += 1
        
        #============================
        # Parse number of bonds
        #============================
        m = re.match(r"\[interactions\]=(\d+)", self.lines[idx])
        if not m:
            raise ValueError(f"Expected [interactions]=<p>, got: {self.lines[idx]!r}")
        self.n_bonds = int(m.group(1))
        idx += 1
        
        #============================
        # Parse bonds
        #============================
        self.bonds = []
        pattern_wo_coupling = r"\(\s*(\d+)\s*,\s*(\d+)\s*\)" # (i, j)
        pattern = r"\s*(\S+)\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)" # J (i, j)
        for _ in range(self.n_bonds):
            bond = {}
            m = re.match(pattern, self.lines[idx])
            if not m:
                m = re.match(pattern_wo_coupling, self.lines[idx])
                if not m:
                    raise ValueError(f"Expected '<J> (<i>, <j>)' or '(<i>, <j>)'', got: {self.lines[idx]!r}")
                bond['coupling'] = ''
                site_i = int(m.group(1))
                site_j = int(m.group(2))
                bond['sites'] = np.array([min(site_i, site_j), max(site_i, site_j)], dtype=int)
            else:
                bond['coupling'] = m.group(1)
                site_i = int(m.group(2))
                site_j = int(m.group(3))
                bond['sites'] = np.array([min(site_i, site_j), max(site_i, site_j)], dtype=int)
            self.bonds.append(bond)
            idx += 1
        
        # safety check - verify that we reached the end of self.lines
        if not idx == len(self.lines):
            raise RuntimeError(f'Problem: parsed {self.n_bonds} bonds, but we did not reach the end of the lattice file ({len(self.lines)} lines).')
        
        # safety check - verify that all bonds are unique
        bonds_l = []
        for bond in self.bonds:
            bonds_l.append(bond['sites'])
        bonds_l_unique = np.unique(bonds_l, axis=0)
        if not len(bonds_l_unique) == len(bonds_l):
            raise RuntimeError(f'Problem: parsed {self.n_bonds} bonds, but there are duplicates.')
        
        #============================
        # Extract couplings
        #============================
        self.couplings = []
        for bond in self.bonds:
            if not bond['coupling'] in self.couplings:
                self.couplings.append(bond['coupling'])
        
        return
    
    
    @staticmethod
    def _get_distance(site1, site2):
        """
        Compute the distance bewteen two sites
    
        Parameters
        ----------
        site1 : numpy array
            coordinates of first site
        site2 : numpy array
            coordinates of second site
    
        Returns
        -------
        d : float
            distance bewteen site1 and site2
        
        """
        d = 0.0
        for i in range(len(site1)):
            d += (site1[i] - site2[i])**2
        d = math.sqrt(d)
        return d
    
    
    def _get_ghosts(self):
        """
        Extract ghost sites for bonds which wrap around the torus.

        """
        
        ghosts = []
        ind_ghosts = []
        
        for bond in self.bonds:
            site1_index = bond['sites'][0]
            site2_index = bond['sites'][1]
            
            site1_coords = self.sites[site1_index][0:2]
            site2_coords = self.sites[site2_index][0:2]
            
            d = Lattice._get_distance(site1_coords, site2_coords)
            
            if abs(d - self.dc[bond['coupling']]) < 1.0e-10:
                bond['site1_mapped_coords'] = site1_coords
                bond['site2_mapped_coords'] = site2_coords
            else:
                # this bond wraps around the torus
                vecs = [self.t1, self.t2,
                        self.t1 + self.t2,
                        self.t1 - self.t2]
                found = False
                # we attempt to find the ghost position of site1
                for j in range(len(vecs)):
                    site1p_coords = site1_coords + vecs[j]
                    d = Lattice._get_distance(site1p_coords, site2_coords)
                    if abs(d - self.dc[bond['coupling']]) < 1.0e-10:
                        self.minX = min(self.minX, site1p_coords[0])
                        self.maxX = max(self.maxX, site1p_coords[0])
                        self.minY = min(self.minY, site1p_coords[1])
                        self.maxY = max(self.maxY, site1p_coords[1])
                        # add ghost point
                        ghosts.append(site1p_coords)
                        ind_ghosts.append(site1_index)
                        bond['site1_mapped_coords'] = site1p_coords
                        bond['site2_mapped_coords'] = site2_coords
                        found = True
                        break
                if found == False:
                    # we attempt to find the ghost position of site2
                    for j in range(len(vecs)):
                        site2p_coords = site2_coords + vecs[j]
                        d = Lattice._get_distance(site1_coords, site2p_coords)
                        if abs(d - self.dc[bond['coupling']]) < 1.0e-10:
                            self.minX = min(self.minX, site2p_coords[0])
                            self.maxX = max(self.maxX, site2p_coords[0])
                            self.minY = min(self.minY, site2p_coords[1])
                            self.maxY = max(self.maxY, site2p_coords[1])
                            # add ghost point
                            ghosts.append(site2p_coords)
                            ind_ghosts.append(site2_index)
                            bond['site1_mapped_coords'] = site1_coords
                            bond['site2_mapped_coords'] = site2p_coords
                            found = True
                            break
                if found == False:
                    print('site1_index = ', site1_index)
                    print('site2_index = ', site2_index)
                    print('site1_pos = ', site1_coords)
                    print('site2_pos = ', site2_coords)
                    print('coupling = ', bond['coupling'])
                    print('target distance = ', self.dc[bond['coupling']])
                    print('base distance = ', Lattice._get_distance(site1_coords, site2_coords))
                    raise RuntimeError('Problem: did not find image point')
        
        ghosts = np.array(ghosts)
        ghosts, ind = np.unique(ghosts, return_index=True, axis=0)
        
        for i in range(len(ghosts)):
            ghost_dict = {}
            ghost_dict['index'] = ind_ghosts[ind[i]]
            ghost_dict['coords'] = ghosts[i]
            self.ghosts.append(ghost_dict)
        
        return
    
    
    def plot(self, **kwargs):
        """
        Plot the lattice
    
        Parameters
        ----------
        show_numbering : bool [optional][default: True]
            If True, plot the numbering of the sites
        plot_bonds : bool [optional][default: True]
            If True, plot the bonds
        plot_simu_torus : bool [optional][default: True]
            If True, plot the simulation torus
        couplings : list [optional]
            list of couplings to plot
        
        Returns
        -------
        fig, ax : matplotlib handles
        
        """
        
        zorder_simu_torus = 1
        zorder_bonds      = 2
        zorder_sites      = 5
        fs_sites = 12
        cols = ['b', 'r', 'm']
        
        origin = np.array([0.0, 0.0])
        
        ax = kwargs.get('ax', None)
        fig = kwargs.get('fig', None)
        show_numbering = kwargs.get('show_numbering', True)
        plot_bonds = kwargs.get('plot_bonds', True)
        plot_simu_torus = kwargs.get('plot_simu_torus', True)
        couplings_to_plot = kwargs.get('couplings', self.couplings)
        
        cols_by_coupling = {}
        cpt = int(0)
        for coupling in self.couplings:
            cols_by_coupling[coupling] = cols[cpt]
            cpt += 1
        
        if (fig == None) & (ax == None):
            fig, ax = plt.subplots(1, 1, figsize=(7, 8), constrained_layout=True)
        else:
            if ((fig == None) & (not ax == None)) | ((not fig == None) & (ax == None)):
                raise RuntimeError('Problem with figure handles. Provide fig and ax.')
        
        #========================================
        # PLOT SIMULATION TORUS
        #========================================
        
        if plot_simu_torus:
            self.minX = min(self.minX, self.t1[0])
            self.minX = min(self.minX, self.t2[0])
            self.maxX = max(self.maxX, self.t1[0])
            self.maxX = max(self.maxX, self.t2[0])
            self.minY = min(self.minY, self.t1[1])
            self.minY = min(self.minY, self.t2[1])
            self.maxY = max(self.maxY, self.t1[1])
            self.maxY = max(self.maxY, self.t2[1])
            ax.quiver(*origin, *self.t1, 
                      angles='xy', scale_units='xy', 
                      scale=1, color='g',
                      zorder=zorder_simu_torus)
            ax.quiver(*origin, *self.t2, 
                      angles='xy', scale_units='xy', 
                      scale=1, color='g',
                      zorder=zorder_simu_torus)
        
        #========================================
        # PLOT LATTICE SITES
        #========================================
        
        ax.scatter(self.sites_x, self.sites_y,
                   marker='o', 
                   s=32,
                   color='k',
                   zorder=zorder_sites,
                   label='_nolegend_')
        
        #========================================
        # PLOT GHOST SITES
        #========================================
        
        sites_ghost_x = np.array([ghost['coords'][0] for ghost in self.ghosts])
        sites_ghost_y = np.array([ghost['coords'][1] for ghost in self.ghosts])
        
        ax.scatter(sites_ghost_x, sites_ghost_y,
                   marker='s', 
                   s=32,
                   color='grey',
                   zorder=zorder_sites,
                   label='_nolegend_')
        
        #========================================
        # PLOT BONDS
        #========================================
        
        if plot_bonds:
            coup = []
            for bond in self.bonds:
                site1_mapped_coords = bond['site1_mapped_coords']
                site2_mapped_coords = bond['site2_mapped_coords']
                coupling = bond['coupling']
                label = '_nolegend_'
                if coupling in couplings_to_plot:
                    if not coupling == '':
                        if not coupling in coup:
                            coup.append(coupling)
                            label = coupling
                    ax.plot(np.array([site1_mapped_coords[0], site2_mapped_coords[0]]),
                            np.array([site1_mapped_coords[1], site2_mapped_coords[1]]),
                            color=cols_by_coupling[coupling],
                            linestyle='--',
                            zorder=zorder_bonds,
                            label=label)
        
        #========================================
        # SHOW NUMBERING
        #========================================
        
        if show_numbering:
            for site_index in self.sites.keys():
                x = self.sites[site_index][0]
                y = self.sites[site_index][1]
                ax.annotate(str(site_index), xy=(x, y),
                            xytext=(5, 5),
                            textcoords='offset points',
                            fontsize=fs_sites)
            
            for ghost in self.ghosts:
                x = ghost['coords'][0]
                y = ghost['coords'][1]
                ghost_index = ghost['index']
                ax.annotate(str(ghost_index), xy=(x, y),
                            xytext=(5, 5),
                            textcoords='offset points',
                            fontsize=fs_sites)
            
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_xlim(left=self.minX-0.4, right=self.maxX+0.4)
        ax.set_ylim(bottom=self.minY-0.4, top=self.maxY+0.4)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_aspect('equal')
        ax.legend(loc='upper right')
        
        return fig, ax
