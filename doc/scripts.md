Sample Scripts    {#scripts}
==============

Liquid Water
----------------------------------------------------------------------------------------------------

The following script can be used to build a simulation box with 500 water molecules and density of
1 g/cm³. In the end, it creates a LAMMPS data file that can be used to run an MD simulation. By
changing the value defined in the first command, one can choose among various 3-site water models.

### File: water.playmol

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define		model as TIP3P
include		$model.playmol
atom_type	HW       $LJ_H
atom_type	OW       $LJ_O
mass		HW       1.008
mass		OW       15.9994
bond_type	HW OW    $bond_OH
angle_type	HW OW HW $angle_HOH
atom		H1 HW    $charge_H
atom		O  OW    $charge_O
atom		H2 HW    $charge_H
bond		O  H1 H2
box    		density 0.602214 # Daltons per cubic angstrom
packmol		diameter 3.0 retry 0.95 pack mol(O) 500 action execute
write		lammps water_$model.lmp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A file that contains model parameters must conform to the LAMMPS model equations. In the case of
harmonic potentials for bond stretching and angle bending, these equations are, respectively:

\f[
  V_{bond} = k_b \left(r - r_0\right)^2 \qquad \text{and} \qquad
  V_{angle} = k_a \left(\theta - \theta_0\right)^2
\f]

The intermolecular potential is:

\f[
  V_{ij} = 4 \epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} - 
                                  \left(\frac{\sigma_{ij}}{r}\right)^6 
                           \right] + C\frac{q_i q_j}{r}
\f]
where \f$\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}\f$ and \f$\sigma_{ij} = (\sigma_i +
\sigma_j)/2\f$. In the files shown below, the units of measurement of the parameters are in
conformance with the LAMMPS unit style `real`.


### File: TIP3P.playmol

This file contains the TIP3P \cite Jorgensen_1983 model parameters:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define LJ_H      as  0.0000 0.0000
define LJ_O      as  0.1520 3.1507
define bond_OH   as  553.00 0.9572
define angle_HOH as  100.00 104.52
define charge_H  as  0.417
define charge_O  as -0.834
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file contains the SPC/Fw \cite Wu_2006 model parameters:

### File: SPC_Fw.playmol

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atom_type   HW        0.0          0.0
atom_type   OW        0.1554253    3.165492
bond_type   HW OW     {1059.162/2} 1.012
angle_type  HW OW HW  {75.90/2}    113.24
charge      OW*      -0.820
charge      HW*       0.410
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

