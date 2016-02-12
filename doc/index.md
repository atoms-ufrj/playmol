Introduction      {#mainpage}
============

Playmol is a(nother) software for building molecular models.

Its main distinguishing features are:

* Molecules are created with simple scripts consisting of a small set of [commands](commands.html).
* Molecular topology arises naturally when atoms are connected (automatic detection of angles and
dihedrals).
* Multiple copies of a molecule are automatically created when new coordinates are defined for their
atoms.
* Integration with [Packmol] \cite Martinez_2009 provides a way of creating complex molecular
systems.
* Generation of [LAMMPS] \cite Plimpton_1995 configuration files provides a way of performing
efficient MD simulations.

The sections of this manual are:

* [Installation](install.html): provides instructions for Playmol installation
* [Playmol Basics](basics.html): explains the working logic of Playmol
* [Commands](commands.html): describes all Playmol commands
* [Sample Scripts](scripts.html): provides sample Playmol scripts

@copyright GNU Public License.

@author Charlles R. A. Abreu (abreu@eq.ufrj.br)<br>
Applied Thermodynamics and Molecular Simulation Group ([ATOMS])<br>
Federal University of Rio de Janeiro, Brazil

[Packmol]:	http://www.ime.unicamp.br/~martinez/packmol
[LAMMPS]:	http://lammps.sandia.gov
[ATOMS]:	http://atoms.peq.coppe.ufrj.br
