Sample Scripts    {#scripts}
==============

TIP3P Water
-----------

The following script builds a simulation box with 1 g/cm³ density and 500 water molecules, and creates a LAMMPS data file to be used for an MD simulation with the TIP3P model parameters:

	# Modified TIP3P model: Price and Brooks, J Chem Phys, 121, 10096 (2004).

	atom_type	H	0.000 0.000
	atom_type	O	0.102 3.188

	mass		H	1.008
	mass		O	15.9994

	bond_type	H O  	450.0 0.9572
	angle_type	H O H	55.0 104.52

	atom		H1	H
	atom		O	O
	atom		H2	H

	charge		H?	 0.415
	charge		O	-0.830

	bond		O 	H1	H2

	box    		density 0.602214 # Daltons per cubic angstrom
	packmol		tolerance 3.0 retry 0.95 pack 1 500 action execute
	write		lammps water.lmp
	write		summary

