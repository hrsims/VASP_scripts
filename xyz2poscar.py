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

# This script is intended to convert an xyz file (particuarly, one in the format saved
# by VESTA) to a VASP5-readable POSCAR. The main use case I envision is creating some
# complicated supercell in VESTA that cannot be readily exported directly (perhaps 
# you imported several phases).

# Lattice parameters for input xyz
a = 50.765
b = 50.765
c = 22.8
# Lattice parameters for new file
a1 = a
b1 = b
c1 = c
# Direct or Cartesian?
dc_string = 'Direct'

# Name of xyz file to convert
infile = 'temp.xyz'
# Name of POSCAR file to write
outfile = 'POSCAR.RT13.vasp'

# End of user edittable portion (...mostly)

# parse xyz file
with open(infile,'r') as f:
    types = []
    ntypes = []
    pos = []
    # This would just be a single line, but I wanted to keep track of the number
    # of atoms of each type
    # Note that this will screw up if your xyz file does not place all of the atoms of
    # the same type together
    for i,line in enumerate(f):
        if i==0:
            natoms = int(line.strip())
        elif i==1:
            label = line
        else:
            temp = line.strip().split()
            t_temp = temp[0]
            if t_temp not in types:
                types.append(t_temp)
                ntypes.append(0)
            ind = types.index(t_temp)
            ntypes[ind] += 1
            pos.append([float(temp[1]),float(temp[2]),float(temp[3])])

# write POSCAR
# this part is not currently general, but xyz files do not contain lattice info
# one could imagine reading another POSCAR for this purpose
# The biggest issue is that it currently assumes no less than orthorhombic
# lattice vectors. That's easy to change, though.

with open(outfile,'w') as f:
    f.write(label)
    f.write("1.0\n")
    f.write("{0:.10f}  {1:.10f}  {2:.10f}\n".format(a,0.,0.))
    f.write("{0:.10f}  {1:.10f}  {2:.10f}\n".format(0.,b,0.))
    f.write("{0:.10f}  {1:.10f}  {2:.10f}\n".format(0.,0.,c))
    for t in types:
        f.write(' '+t+'  ')
    f.write('\n')
    for n in ntypes:
        f.write('{0:d}  '.format(n))
    f.write('\n')
    f.write(dc_string+'\n')
    # If your two crystals don't perfectly match, you might want to include checks
    # here to catch stray atoms outside the unit cell. I dunno!
    if dc_string=='Direct':
        for p in pos:
            f.write('{0:.10f}  {1:.10f}  {2:.10f}\n'.format(p[0]/a1,p[1]/b1,p[2]/c1))
    else:
        for p in pos:
            f.write('{0:.10f}  {1:.10f}  {2:.10f}\n'.format(p[0],p[1],p[2]))

