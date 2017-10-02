# Sucrose Example

----------------------------------------------------------------------------------------------------

![sucrose](sucrose.png)
> Visualization created using [VMD](www.ks.uiuc.edu/Research/vmd/)

----------------------------------------------------------------------------------------------------
Requirements
----------------------------------------------------------------------------------------------------

1. [AMBER Tools](http://ambermd.org/#AmberTools) installed
2. Environmental variable AMBERHOME containing the path to AMBER Tools' directory

----------------------------------------------------------------------------------------------------
Preparation Steps (Done!)
----------------------------------------------------------------------------------------------------

Create file `GLYCAM_06j.playmol` from AMBER's `GLYCAM_06j.dat` using playmoltools:

    playmoltools -f amber -i $AMBERHOME/dat/leap/parm/GLYCAM_06j.dat -o GLYCAM_06j.playmol

Create file `sucrose.pdb` using Glycam-Web's [Carbohydrate Builder](http://glycam.org/cb). For this,
use the condensed code `DFrufb2-1DGlcpa` to generate a PDB file, download it, and execute:

    mv 1.pdb sucrose.pdb

Create file `sucrose.playmol` from `sucrose.pdb` using playmoltools and GLYCAM's prep file:

    playmoltools -f pdb -p $AMBERHOME/dat/leap/prep/GLYCAM_06j-1.prep -i sucrose.pdb -o sucrose.playmol

----------------------------------------------------------------------------------------------------
Execution Step
----------------------------------------------------------------------------------------------------

Run Playmol to create a system containing one sucrose molecule solvated in water:

    playmol sucrose_in_water.playmol

