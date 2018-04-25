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

def read_poscar_header():
    """Grab information from POSCAR about ionic species"""
    infile = "POSCAR"
    try:
        with open(infile,'r') as f:
            # We're not worried about the lattice, so skip this info
            for i in range(5):
                f.readline()
            labels = f.readline().strip().split()
            num = [np.float(n) for n in f.readline().strip().split()]
    except:
        print("Could not open file POSCAR")
        exit(1)
    return labels,num

def write_proj_files(energies,proj_array,proj_array_tot,nspecies,nspin,nk,nbands,which_orb):
    with open("proj_out.dat","w") as f:
        # first write info to be read by proplot
        f.write('{0:d} {1:d} {2:d} {3:d}\n'.format(nspecies,nspin,nk,nbands))
        for ispin in range(nspin):
            for ik in range(nk):
                for ib in range(nbands):
                    f.write("{0:.8f}  ".format(energies[ispin,ik,ib]))
            f.write("\n")
        f.write("\n")
        for isp in range(nspecies):
            for ik in range(nk):
                for ib in range(nbands):
                    f.write("{0:.8f}  ".format(proj_array_total[isp,ispin,ik,ib,which_orb[isp]]))
                f.write("\n")
            f.write("\n")

    with open("proj_out_total.dat","w") as f:
        # first write info to be read by proplot
        f.write('{0:d} {1:d} {2:d} {3:d}\n'.format(nspecies,nspin,nk,nbands))
        for ispin in range(nspin):
            for ik in range(nk):
                for ib in range(nbands):
                    f.write("{0:.8f}  ".format(energies[ispin,ik,ib]))
            f.write("\n")
        f.write("\n")
        for isp in range(nspecies):
            for ik in range(nk):
                for ib in range(nbands):
                    f.write("{0:.8f}  ".format(proj_array_total[isp,ispin,ik,ib]))
                f.write("\n")
            f.write("\n")

if __name__ == "__main__":
    """Read in PROCAR, read and store the orbital character of the bands from one
    or more ions, and write the data in a reasonably useful form, particularly if
    you are using my plotting script."""
    infile = "PROCAR"
    parser = argparse.ArgumentParser(description='Read a VASP PROCAR file and write'+
                                     'useful collections of projections.')
    parser.add_argument('--nspin',type=int,help='number of spin channels (default 1)')
    parser.add_argument('--orbs',nargs='+',help='list of which orbital to print per species (s,p,d) (default d)')

    args = parser.parse_args()

    nspin = args.nspin
    which_orb_c = args.orbs
    if nspin is None: nspin = 1

    types, ntypes = read_poscar_header()
    nspecies = len(types)
    if which_orb_c is None: which_orb_c = ['d']*nspecies
    orbdict = {'s':0,'p':1,'d':2}
    which_orb = []
    for o in which_orb_c:
        which_orb.append(orbdict[o])

    try:
        with open(infile,'r') as f:
            f.readline()
            nums = f.readline().strip().split()
            # number of k-points, bands, and ions
            nk = int(nums[3])
            nbands = int(nums[7])
            nions = int(nums[11])
            # store the band energies and the orbital projections (s,p,d)
            energies = np.zeros((nspin,nk,nbands),dtype=float)
            proj_array = np.zeros((nspecies,nspin,nk,nbands,3),dtype=float)
            # maybe you also want the summed contributions from each species
            proj_array_tot = np.zeros((nspecies,nspin,nk,nbands),dtype=float)
            # Let's get started.
            # SO MANY for loops, but I'm reading a file...
            for ispin in range(nspin):
                if (ispin==1):
                    f.readline() # blank
                for ik in range(nk):
                    f.readline() # blank
                    f.readline() # which kpt and what weight?
                    for ib in range(nbands):
                        f.readline() # blank
                        energies[ispin,ik,ib] = float(f.readline().strip().split()[4])
                        f.readline() # blank
                        f.readline() # labels
                        for isp in range(nspecies):
                            for i in range(ntypes[ispin]):
                                temp = f.readline().strip().split()
                                proj_array[isp,ispin,ik,ib,0] = float(temp[1])
                                proj_array[isp,ispin,ik,ib,1] = sum([float(x) for x in temp[2:5]])
                                proj_array[isp,ispin,ik,ib,2] = sum([float(x) for x in temp[5:10]])
                                # in case you also want all contributions from a species
                                proj_array_tot[isp,ispin,ik,ib] += float(temp[10])
                            temp = f.readline().strip().split()
                            # Normalize by the total weight, which is generally a bit less than 1
                            proj_array[isp,ispin,ik,ib,:] /= float(temp[10]) # total projections
                            proj_array_tot[isp,ispin,ik,ib] /= float(temp[10]) # total projections
                    f.readline()
    except:
        print("Cannot open file "+infile)
        exit(1)
    write_proj_files(energies,proj_array,proj_array_tot,nspecies,nspin,nk,nbands,which_orb)

