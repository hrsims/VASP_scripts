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

def setup_kline(lens,length):
    """Completely non-general, and the labels must be added by hand :("""
    x=[0,lens[0],sum(lens[:2]),sum(lens[:3]),sum(lens)]
    labels=['-M','-K',r'$\Gamma$','K','M']
    return x,labels

def path_length(ki,kf,a):
    pref = 2*np.pi/a
    return pref*np.linalg.norm(kf-ki)

def set_kpoints():
    kpts = []
    try:
        with open("KPOINTS_bands",'r') as f:
            f.readline() # comment
            ndiv = int(f.readline().strip())
            f.readline() # had better be Line
            RorC = f.readline().strip()[0] # Reciprocal or Cartesian?
            for line in f:
                if line.strip() != '':
                    kpts.append([float(x) for x in line.strip().split()])
    except:
        print("Cannot open file KPOINTS_bands")
        exit(1)
    kpts = np.array(kpts)
    len_array = []
    x_array = []
    for i,k in enumerate(kpts):
        len_array.append(path_length(kpts[i-1],kpts[i],a))
        x_array.append(np.linspace(sum(len_array[:-1]),len_array[-1],ndiv,endpoint=True))
    fulllen = sum(len_array)
    xx = np.concatenate(x_array)
    x,labels = setup_kline(len_array,fulllen)
    return xx,x,labels

if __name__ == "__main__":
    # First open the file just to get some info
    with open('proj_out.dat','r') as f:
        nspecies,nspin,nk,nbands = [int(n) for n in f.readline().strip().split()]

    parser = argparse.ArgumentParser(description='Read projections from proj_out.dat'+
                                     'and plot them.')
    parser.add_argument('--ef',metavar='Ef',type=int,help='Fermi energy in eV (no default)')
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
        for ispin in range(nspin):
            for ik in range(nk):
                bandline = f.readline().strip().split()
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

    # This should not have to be entered by hand!
    a = 5.431
    kpts,x,labels = set_kpoints()

    fig, ax = plt.subplots(1,nspin)
    for isp in range(nspecies):
        for ispin in range(nspin):
            for i in range(nbands):
                col = proj_array[isp,ispin,:,i]
                ax[ispin].scatter(range(nk),np.array(energies[0,:,i])-ef,\
                      c=col,cmap=plt.cm.Reds,s=10*col,edgecolor=None)
    ax[:].set_ylim(emin,emax)
    ax[:].set_xlim(0,nk-1)
    ax[:].axhline(0,c='k',ls='--')
    ax[:].set_ylabel(r"$E - E_F$ (eV)")
    ax[:].xaxis.set_ticks(x)
    ax[:].xaxis.set_ticklabels(labels)
    for xx in x:
        ax[:].axvline(xx,c='k')
    plt.savefig(outfile+".png",dpi=300)
    plt.show()
