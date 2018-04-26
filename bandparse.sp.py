#! /usr/bin/env python

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

def path_length(ki,kf,a,RorC):
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
                    kpts.append([float(x) for x in line.strip().split()[:3]]) # allow for comments
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
    with open('EIGENVAL','r') as f:
      for i in range(5):
        f.readline()
      dummy, nk, n = [int(x) for x in f.readline().split()]
      f.readline()
      e = [[] for i in range(n)]
      kx = [[] for i in range(nk)]
      ky = [[] for i in range(nk)]
      kz = [[] for i in range(nk)]
      for i in range(nk):
        kx[i],ky[i],kz[i],wgt = [float(x) for x in f.readline().split()]
        for j in range(n):
          temp = f.readline().split()
          ii, enup, endn = int(temp[0]), float(temp[1]), float(temp[2])
          e[ii-1].append([enup,endn])
        f.readline()

    parser = argparse.ArgumentParser(description='Plot the band structure from a VASP EIGENVAL file.')
    parser.add_argument('file_name',help='file name to save png band plot (default bands)')
    parser.add_argument('--nspin',help='number of spin channels (default 1)')
    parser.add_argument('--ef',metavar="Ef",type=float,help='Fermi Energy in eV (no default)')
    parser.add_argument('--emin',metavar='Emin',type=float,help='minimum energy in plot window (default -5)')
    parser.add_argument('--emax',metavar='Emax',type=float,help='maximum energy in plot window (default 5)')

    args = parser.parse_args()

    file_name = args.file_name
    nspin = args.nspin
    ef = args.ef
    emin = args.emin
    emax = args.emax
    if file_name is None:
        bfile = 'bands.png'
    else:
        bfile = file_name+'png'
    if nspin is None: nspin = 1
    if emin is None: emin = -5
    if emax is None: emax = 5

    e = np.array(e)

    fig, ax = plt.subplots(1,nspin)
    xx,x,labels = setup_kline(x,labels)
    colors = {0:'k',1:'r'}
    for ispin in range(nspin):
        for i in range(n):
            ax[ispin].scatter(xx,np.array(e[i,:,ispin])-ef,c=colors[ispin])
        ax[ispin].set_ylim(emin,emax)
        ax[ispin].set_xlim(0,max(xx))
        ax[ispin].axhline(0,c='k',ls='--')
        ax[ispin].xaxis.set_ticks(x,labels)
    ax[0].set_ylabel(r"$E - E_F$ (eV)")
    plt.savefig(bfile,dpi=300)
    plt.show()
