Playmol tools
=============

----------------------------------------------------------------------------------------------------
amber2playmol
----------------------------------------------------------------------------------------------------

This tool converts an [AMBER] parameter file into a Playmol script.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
amber2playmol.py < /opt/amber17/dat/leap/parm/gaff.dat > gaff.playmol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

----------------------------------------------------------------------------------------------------
pqr2playmol
----------------------------------------------------------------------------------------------------

This tool converts a PQR file into a Playmol script. The format of a PQR file is similar to the
[PDB file format], but with _charge_ and _radius_ replacing the fields _occupancy_ and _tempFactor_,
respectively. The fields are:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
recordName serial atomName residueName chainID residueNumber x y z charge radius atomType
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PQR files can be created using [Avogadro] by adding the extension _.pqr_ when saving a file after
choosing `Save as...` in the menu entry `file`. [AMBER]'s package [Antechamber] is also able to save
files in the PQR format by using the option `-fo mpdb`.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pqr2playmol.py < ibuprofen.pqr > ibuprofen.playmol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


[AMBER]:		http://ambermd.org
[PDB file format]:	http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
[Avogadro]:		http://avogadro.cc
[Antechamber]:		http://ambermd.org/antechamber/ac.html
