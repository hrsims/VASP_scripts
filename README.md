# VASP_scripts
Some useful (mostly Python) scripts I use for postprocessing VASP. Generally not overlapping with the very nice VTST ones.

bandparse.sp.py
A script to read the EIGENVAL file and produce a band plot. Sp means that it supports spin-polarized bands.

proparse_multi.py
A script to read the PROCAR file and produce a data file containing both the band structure and the projection of the wave function on atomic wave functions at the atomic sites.

proplot_multi.py
A script to read the output from proparse_multi.py and plot it in a "fat bands" manner, i.e. such that the color and size of the band is related to the orbital character at that point.

xyz2poscar.py
A script to convert an xyz file derived from a known structure (perhaps created in a way that makes it difficult or impossible to save directly as a POSCAR) to a POSCAR format file.

POSCAR, KPOINTS_bands, EIGENVAL, PROCAR
Sample VASP output for Si
