Playmol
=======

Playmol is a(nother) software for building molecular models.

Its main distinguishing features are:

* Molecules are created with simple scripts consisting of a small set of commands.
* Molecular topology arises naturally when atoms are connected (automatic detection of angles and
dihedrals).
* Multiple copies of a molecule are automatically created when new coordinates are defined for their
atoms.
* Integration with [Packmol] provides a way of creating complex molecular systems.
* Generation of [LAMMPS] configuration files provides a way of performing efficient MD simulations.

Copyright: GNU Public License.

Author: [Charlles R. A. Abreu](mailto:abreu@eq.ufrj.br)<br>
Applied Thermodynamics and Molecular Simulation Group ([ATOMS])<br>
Federal University of Rio de Janeiro, Brazil

[Packmol]:	http://www.ime.unicamp.br/~martinez/packmol
[LAMMPS]:	http://lammps.sandia.gov
[ATOMS]:	http://atoms.peq.coppe.ufrj.br

--------------------------------------------------------------------------------

Online Manual
-------------

The Playmol User Manual is available at:

[http://atoms-ufrj.github.io/playmol/docs/html](http://atoms-ufrj.github.io/playmol/docs/html)

--------------------------------------------------------------------------------

Forum
-----

A forum for Questions & Answers about Playmol is available at:

[https://groups.google.com/forum/#!forum/playmol-users](https://groups.google.com/forum/#!forum/playmol-users)

--------------------------------------------------------------------------------

Installation
------------

Playmol is distributed as a git repository. To download it, just run:

    git clone https://github.com/atoms-ufrj/playmol

To compile the source code and install Playmol in your system, you can do:

    cd playmol
    make
    sudo make install

To update Playmol, enter the playmol directory and execute the following commands (including
recompilation and reinstallation):

    git pull
    make
    sudo make install

--------------------------------------------------------------------------------

User's Manual
-------------

The Playmol User's Manual is available online [here](http://atoms.peq.coppe.ufrj.br/playmol). You
can also generate a local version if you have [Doxygen](http://www.doxygen.org) (version 1.8 or
later) installed in your system. If you do not have Doxygen, you can download and install it by:

    wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.11.src.tar.gz
    tar -zxvf doxygen-1.8.11.src.tar.gz
    cd doxygen-1.8.11/
    sudo apt-get install cmake flex bison
    mkdir build && cd build/
    cmake ../
    make
    sudo make install

In order to generate the local User's Manual, please go to the playmol directory and execute:

    make doc

The manual will be available as a file _playmol/doc/html/index.html_, which you can open using your
favorite web browser.

--------------------------------------------------------------------------------

Using Playmol
-------------

Once Playmol is installed, you can execute a series of input scripts by typing:

    playmol file-1 [file-2 ...]

This will execute the files in sequence as if they were a unique script. To execute the scripts one
at a time, just run playmol multiple times.

Another way of runnig a playmol script is by starting it with the following line and then making it
executable (e.g. via chmod +x):

    #!/usr/local/bin/playmol

--------------------------------------------------------------------------------

List of Playmol Commands
------------------------

Here is a complete list of Playmol commands:

| Command           | Description                                                               |
|:-----------------:|:--------------------------------------------------------------------------|
| **define**        | defines a variable for further substitution                               |
| **for/next**      | executes commands repeatedly while changing the value of a variable       |
| **if/then/else**  | executes commands conditionally or selects between two command sequences  |
| **atom_type**     | creates an atom type with given name and parameters                       |
| **mass**          | specifies the mass of atoms of a given type                               |
| **diameter**      | specifies the diameter of atoms of a given type                           |
| **bond_type**     | defines parameters for bonds between atoms of two given types             |
| **angle_type**    | defines parameters for angles involving atoms of three given types        |
| **dihedral_type** | defines parameters for dihedrals involving atoms of four given types      |
| **improper_type** | defines parameters for impropers involving atoms of four given types      |
| **models**        | defines whether the first attribute of each type is a model name          |
| **atom**          | creates an atom with given name and type                                  |
| **charge**        | specifies the charge of a given atom                                      |
| **bond**          | creates chemical bonds and automatically detects angles and dihedrals     |
| **improper**      | creates an improper involving four given atoms or search for impropers    |
| **body**          | selects a set of atoms to be marked as a single body                      |
| **extra**         | creates an extra bond, angle, or dihedral involving given atoms           |
| **link**          | links two atoms (and fuses their molecules) without actually bonding them |
| **unlink**        | removes an existing link (and splits the corresponding molecule)          |
| **build**         | guesses atom positions from provided geometric information                |
| **prefix/suffix** | defines default prefixes or suffixes for atom types and atoms             |
| **box**           | defines the properties of a simulation box                                |
| **velocity**      | defines parameters for generating Maxwell-Boltzmann atomic velocities     |
| **align**         | aligns the principal axes of a molecule to the Cartesian axes             |
| **packmol**       | executes Packmol in order to create a packed molecular system             |
| **write**         | saves system info in different file formats (including LAMMPS data files) |
| **include**       | includes commands from another script                                     |
| **reset**         | resets a list of entities together with its dependent lists               |
| **shell**         | executes an external shell command                                        |
| **quit**          | interrupts the execution of a Playmol script                              |

The syntax and behavior of each command is described in the Playmol documentation.

--------------------------------------------------------------------------------

Input Script Examples
-------------------------

Some input script examples are available in the playmol/examples directory.

