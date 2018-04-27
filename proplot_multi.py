# Copyright (C) Hunter Sims 2018
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import numpy as np
import matplotlib.pyplot as plt

def setup_k_label(lens):
    """Sadly the labels must still be added by hand :("""
    x = [0,]
    for l in lens:
        x.append(l+x[-1])
    labels=['L',r'$\Gamma$','X',r'$\Gamma$']
    return x,labels

def path_length(ki,kf,a,b,c):
    pref = 2*np.pi/np.array([a,b,c])
    return np.linalg.norm(pref*(kf-ki))

def set_kpoints(a,b,c):
    kpts = []
    try:
        with open("KPOINTS_bands",'r') as f:
            f.readline() # comment
            ndiv = int(f.readline().strip())
            f.readline() # had better be Line
            RorC = f.readline().strip()[0] # Reciprocal or Cartesian?, not currently used
            for line in f:
                if line.strip() != '':
                    kpts.append([float(x) for x in line.strip().split()])
    except:
        print("Cannot open file KPOINTS_bands")
        exit(1)
    # list of special points
    kpts = np.array(kpts)
    len_array = []
    x_array = []
    for i,k in enumerate(kpts):
        # only look at the end of the paths
        if i%2==0:
            continue
        else:
            len_array.append(path_length(kpts[i-1],kpts[i],a,b,c))
        x_array.append(np.linspace(sum(len_array[:-1]),sum(len_array[:-1])+len_array[-1],\
                                   ndiv,endpoint=True))
    fulllen = sum(len_array)
    xx = np.concatenate(x_array)
    x,labels = setup_k_label(len_array)
    return xx,x,labels

def read_latt_const():
    """Grab information from POSCAR about ionic species"""
    infile = "POSCAR"
    try:
        with open(infile,'r') as f:
            f.readline() # comment
            a0 = float(f.readline().strip())
            avec = []
            for i in range(3):
                avec.append([float(x) for x in f.readline().strip().split()])
            if a0==1.0:
                a = np.linalg.norm(avec[0])
                b = np.linalg.norm(avec[1])
                c = np.linalg.norm(avec[2])
            else:
                a = np.linalg.norm(avec[0])*a0
                b = np.linalg.norm(avec[1])*a0
                c = np.linalg.norm(avec[2])*a0
    except:
        print("Could not open file POSCAR")
        exit(1)
    return a,b,c

if __name__ == "__main__":
    # First open the file just to get some info
    with open('proj_out.dat','r') as f:
        nspecies,nspin,nk,nbands = [int(n) for n in f.readline().strip().split()]

    parser = argparse.ArgumentParser(description='Read projections from proj_out.dat'+
                                     'and plot them.')
    parser.add_argument('--ef',metavar='Ef',type=float,help='Fermi energy in eV (no default)')
    parser.add_argument('--emin',metavar='Emin',type=float,help='minimum energy in plot window (default -5)')
    parser.add_argument('--emax',metavar='Emax',type=float,help='maximum energy in plot window (default 5)')
    parser.add_argument('outfile',help='name of (png) file in which to save bands (default bands.png)')

    args = parser.parse_args()

    ef = args.ef
    emin = args.emin
    emax = args.emax
    outfile = args.outfile

    if emin is None: emin = -5
    if emax is None: emax = 5
    if outfile is None: outfile = 'bands'

    energies = np.zeros((nspin,nk,nbands),dtype=float)
    proj_array = np.zeros((nspecies,nspin,nk,nbands),dtype=float)
    with open("proj_out.dat",'r') as f:
        f.readline() # variables
        for ispin in range(nspin):
            for ik in range(nk):
                bandline = f.readline().strip().split()
                if bandline == ['']: continue
                for ib in range(nbands):
                    energies[ispin,ik,ib] = float(bandline[ib])
            f.readline()
        for isp in range(nspecies):
            for ispin in range(nspin):
                for ik in range(nk):
                    bandline = f.readline().strip().split()
                    for ib in range(nbands):
                        proj_array[isp,ispin,ik,ib] = float(bandline[ib])
                f.readline()

    a,b,c = read_latt_const()
    kpts,x,labels = set_kpoints(a,b,c)

    fig, ax = plt.subplots(1,nspin)

    # For more than 5 species, you're on your own
    cmaps = {0:plt.cm.Reds,1:plt.cm.Blues,2:plt.cm.Oranges,3:plt.cm.Purples,4:plt.cm.Greys}
    if hasattr(ax,'__iter__'):
        for ispin in range(nspin):
            for isp in range(nspecies):
                for i in range(nbands):
                    col = proj_array[isp,ispin,:,i]
                    ax[ispin].scatter(kpts,np.array(energies[0,:,i])-ef,\
                                      c=col,cmap=cmaps[isp],s=10*col,edgecolor=None)
            ax[ispin].set_ylim(emin,emax)
            ax[ispin].set_xlim(0,max(kpts))
            ax[ispin].axhline(0,c='k',ls='--')
            ax[ispin].xaxis.set_ticks(x)
            ax[ispin].xaxis.set_ticklabels(labels)
            for xx in x:
                ax[ispin].axvline(xx,c='k')
        ax[0].set_ylabel(r"$E - E_F$ (eV)")
    else:
        for isp in range(nspecies):
            for i in range(nbands):
                col = proj_array[isp,ispin,:,i]
                ax.scatter(kpts,np.array(energies[0,:,i])-ef,\
                           c=col,cmap=cmaps[isp],s=10*col,edgecolor=None)
        ax.set_ylim(emin,emax)
        ax.set_xlim(0,max(kpts))
        ax.axhline(0,c='k',ls='--')
        ax.xaxis.set_ticks(x)
        ax.xaxis.set_ticklabels(labels)
        for xx in x:
            ax.axvline(xx,c='k')
        ax.set_ylabel(r"$E - E_F$ (eV)")

    plt.savefig(outfile+".png",dpi=300)
    plt.show()
