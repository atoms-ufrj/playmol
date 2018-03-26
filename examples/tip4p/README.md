# TIP4P Example

----------------------------------------------------------------------------------------------------

[![TIP4P](http://www.sklogwiki.org/SklogWiki/images/thumb/a/a5/Four_site_water_model.png/300px-Four_site_water_model.png)](http://www.sklogwiki.org/SklogWiki/index.php/TIP4P_model_of_water)

----------------------------------------------------------------------------------------------------
Description
----------------------------------------------------------------------------------------------------

This example uses Playmol, OpenMM, and LAMMPS to simulate a system of rigid, four-site water
molecules.

----------------------------------------------------------------------------------------------------
Provided Files
----------------------------------------------------------------------------------------------------

1. __tip4p.playmol__: Playmol input script
2. __tip4p-2005.params__ and others: Parameter values for various water models. The available models
are:
* TIP4P, [J Chem Phys 79, 926 (1983)](https://doi.org/10.1063/1.445869)
* TIP4P/2005, [J Chem Phys, 123, 234505 (2005)](http://dx.doi.org/10.1063/1.2121687)
* TIP4P-Ew, [J Chem Phys 120, 9665 (2004)](https://doi.org/10.1063/1.1683075)
* TIP4P/Ice, [J Chem Phys 122, 234511 (2005)](https://doi.org/10.1063/1.1931662)
3. __in.lammps__: LAMMPS input script
4. __openmm.py__: OpenMM input script

----------------------------------------------------------------------------------------------------
Requirements
----------------------------------------------------------------------------------------------------

1. [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/linux.html)
2. [OpenMM](http://docs.openmm.org/latest/userguide/application.html#installing-openmm)
3. [LAMMPS](http://lammps.sandia.gov/doc/Section_start.html#making-lammps)
4. [CUDA Toolkit](https://developer.nvidia.com/cuda-90-download-archive) (Optional)

----------------------------------------------------------------------------------------------------
Generation of Initial Configuration Files
----------------------------------------------------------------------------------------------------

Run Playmol to create a system containing a number of water molecules:

    playmol tip4p.playmol

This will create files `tip4p-model.lmp`, `tip4p-model.pdb`, and `tip4p-model.xml`. The first one
is required for running an MD simulation with LAMMPS, while the other two files are required for
running it with OpenMM.

----------------------------------------------------------------------------------------------------
Simulation
----------------------------------------------------------------------------------------------------

Execute LAMMPS to simulate the system:

    mpirun -n 2 $LAMMPSHOME/src/lmp_mpi -in in.lammps

Execute OpenMM to simulate the system:

    python openmm.py

__NOTE__: If CUDA is available, just change line 9 in file `openmm.py` in order to use a GPU for
speeding up the simulation.
