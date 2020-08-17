Sample Scripts    {#scripts}
==============

| Script                | System Description                                                       |
|:---------------------:|:-------------------------------------------------------------------------|
| [Liquid Water]        | liquid water with various force fields                                   |
| [Linear Hydrocarbons] | united-atom, linear hydrocarbons of arbitrary chain length               |

----------------------------------------------------------------------------------------------------
<a name="liquid_water"></a>
Liquid Water
----------------------------------------------------------------------------------------------------

The following script can be used to build a simulation box with 500 water molecules and density of
1 g/cm³. In the end, it creates a LAMMPS data file that can be used to run an MD simulation. By
changing the value defined in the first command, one can choose among various 3-site water models.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define 		model as TIP3P
define		density as 0.998 # g/cm³
define  	N as 500
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
box    		density {0.602214*$density} # Da/Å³
packmol		diameter 3.0 retry 0.95 pack mol(O) $N action execute
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
\sigma_j)/2\f$. In the files shown below, the units of the parameters are in conformance with the
LAMMPS unit style `real`.

### File: TIP3P.playmol

The following file contains the TIP3P \cite Jorgensen_1983 model parameters:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define LJ_H      as  0.0000 0.0000
define LJ_O      as  0.1520 3.1507
define bond_OH   as  553.00 0.9572
define angle_HOH as  100.00 104.52
define charge_H  as  0.417
define charge_O  as -0.834
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### File: SPC_E.playmol

The following file contains the SPC/E \cite Berendsen_1987 model parameters:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define LJ_H      as  0.0000       0.0000
define LJ_O      as  0.1554253    3.165492
define bond_OH   as  {1059.162/2} 1.0        # k from SPC/Fw
define angle_HOH as  {75.90/2}    109.47     # k from SPC/Fw
define charge_H  as  0.4238
define charge_O  as -0.8476
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### File: SPC_Fw.playmol

The following file contains the SPC/Fw \cite Wu_2006 model parameters:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define LJ_H      as  0.0000       0.0000
define LJ_O      as  0.1554253    3.165492
define bond_OH   as  {1059.162/2} 1.012
define angle_HOH as  {75.90/2}    113.24
define charge_H  as  0.410
define charge_O  as -0.820
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



----------------------------------------------------------------------------------------------------
<a name="linear_hydrocarbons"></a>
Linear Hydrocarbons
----------------------------------------------------------------------------------------------------

The following script creates a simulation box for with single-component, united-atom, linear
hydrocarbon molecules. The units are consistent with the LAMMPS unit style `real`.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define N     as 5
define FF    as nerd
define seed  as 3467
define nmols as 500
define rho   as 0.6 # g/cm³
define temp  as 180 # K

include $FF.playmol

mass CH4 16.04242 # Da
mass CH3 15.03450 # Da
mass CH2 14.02658 # Da

if {$N == 1} then
  atom C1 CH4
else
  atom C1 CH3
  for i from 2 to {$N-1}
    atom C$i CH2
    bond C$i C{$i-1}
  next
  atom C$N CH3
  bond C$N C{$N-1}
endif

if {$N > 3} then
  build
    $N
    C1 0 0 0
    C2 C1 $L
    C3 C2 $L C1 $theta
    for i from 4 to $N
      C$i C{$i-1} $L C{$i-2} $theta C{$i-3} $phi
    next
endif

box      density {0.602214*$rho} # Da/Å³
align    mol(C1) x y
packmol  seed $seed retry 0.95 copy mol(C1) $nmols action execute

write    lmp/models C${N}_$FF.lammpsdata
write    summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### File: nerd.playmol

The following file contains the NERD \cite Nath_1998 model parameters:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define          cutoff  as 13.8     # Å
define          L       as 1.54     # Å
define          theta   as 114      # degrees
define          phi     as 180      # degrees
define          kB      as 1.987E-3 # kcal/mol.K

atom_type       CH2     lj/cut  {45.8*$kB} 3.930 $cutoff
if {$N == 2} then
  atom_type     CH3     lj/cut {100.6*$kB} 3.825 $cutoff
else
  if {$N == 3} then
    atom_type   CH3     lj/cut {102.6*$kB} 3.857 $cutoff
  else
    atom_type   CH3     lj/cut {104.0*$kB} 3.910 $cutoff
  endif
endif
atom_type       CH4     lj/cut {148.0*$kB} 3.730 $cutoff # TraPPE parameters

diameter        CH4     3.730 # Å
diameter        CH3     3.910 # Å
diameter        CH2     3.730 # Å

bond_type     CH? CH?         harmonic {96500*$kB/2} $L
angle_type    CH? CH? CH?     harmonic {62500*$kB/2} $theta
dihedral_type CH? CH? CH? CH? opls {355.04*2*$kB} {-68.19*2*$kB} {701.32*2*$kB} 0.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

----------------------------------------------------------------------------------------------------

[Liquid Water]:			#liquid_water
[Linear Hydrocarbons]:		#linear_hydrocarbons
