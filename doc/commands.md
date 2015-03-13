Commands      {#commands}
========

Playmol is designed to execute scripts containing the commands described below:

* [atom_type] - creates an atom type with given name and parameters.
* [mass] - specifies the mass of atoms of a given type.
* [bond_type] - defines parameters for bonds between atoms of two given types.
* [angle_type] - defines parameters for angles involving atoms of three given types.
* [dihedral_type] - defines parameters for dihedrals involving atoms of four given types.
* [improper_type] - defines parameters for impropers involving atoms of four given types.
* [atom] - creates an atom with given name and type.
* [charge] - specifies the charge of a given atom.
* [bond] - creates a bond between two given atoms (angles and dihedrals are automatically detected).
* [improper] - creates an improper involving four given atoms.
* [extra_dihedral] - creates an extra dihedral involving four given atoms.
* [xyz] - defines positions for all atoms of one or more molecules.
* [box] - defines the properties of a simulation box.
* [packmol] - executes Packmol to create a packed molecular system.
* [write] - writes down system info in different file formats (including LAMMPS data files).
* [prefix] - defines default prefixes for atom types and atoms.
* [suffix] - defines default suffixes for atom types and atoms.
* [include] - includes commands from another script.
* [reset] - resets a list of entities together with its dependent lists.
* [shell] - executes an external shell command.
* [quit] - interrupts the execution of a Playmol script.


------------
Introduction
------------

Playmol scripts are text files containing commands whose syntaxes are described in this section. Each command is a sequence of keywords and parameter values separated by spaces and/or tabs. A script can include comments, identified by the comment mark "#". When such mark is found, all trailing characters in the same line (including the comment mark itself) are ignored. A command can span several lines by means of continuation marks "..." or "&". When the actual part of a command line (comments excluded) ends with a continuation mark, the command will continue in the next line (of course, the continuation mark itself is not part of the command).

In the examples of this section, the units employed for physically meaningful values are those corresponding to [LAMMPS real units].

-------------------------------
<a name="atom_type"/> atom_type
-------------------------------

**Syntax**:

	atom_type 	<id> <attribute-list>

* _id_ = a name to be assigned to the atom type being defined
* _attribute-list_ = a list of attributes related to the atom type being defined

**Description**:

This command defines an atom type and its related attributes. At least one atom type must necessarily be defined before any [atom] is created and before any [bond_type], [angle_type], [dihedral_type], or [improper_type] is defined.

The parameter _id_ must be a single string with no comment tags (#) and no wildcard characters (* or ?). If a type [prefix] has been previously activated, then the actual atom type identifier will contain such prefix, followed by the string in _id_. Two atom types cannot have the same identifier.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS pair style associated with the atom type in question.

**Examples**:

	atom_type	O	0.1553 3.166
	atom_type	H	0.0000 0.000

In the example above, atom types _O_ and _H_ are defined with Lennard-Jones parameters _epsilon_ and _sigma_, respectively, corresponding to the SPC water model. This can be used to generate a LAMMPS configuration file associated with the pair style _lj/cut/coul/long_, for instance.

	prefix		type TraPPE-
	atom_type	CH3 lj/cut 0.1947 3.75
	atom_type	CH2 lj/cut 0.0914 3.95

In the example above, united-atom types _TraPPE-CH3_ and _TraPPE-CH2_ are defined with TraPPE Lennard-Jones parameters _epsilon_ and _sigma_, respectively. An explicit pair style specification like _lj/cut_ is necessary in LAMMPS if pair style _hybrid_ or _hybrid/overlay_ is employed.

**See also**:

[atom], [bond_type], [angle_type], [dihedral_type], [improper_type], [prefix]

---------------------
<a name="mass"/> mass
---------------------

**Syntax**:

	mass		<type> <value>

* _type_ = the name of a atom type
* _value_ = mass of any atom with the specified atom type

**Description**:

This command defines the mass of atoms of a given type. Wildcard characters (* or ?) can be used to assign the same mass value to different atom types. The command applies for all atom types, defined either beforehand or afterwards.

An [atom] will not be created if a mass value has not been previously defined to its corresponding atom type.

**Examples**:

	mass		CH2 14.0
	mass		H* 1.008

In the example above, a mass value of _14.0_ is assigned to atom type _CH2_ and a mass value of _1.008_ is assigned to all atoms types whose names start with _H_.

**See also**:

[atom_type], [atom]

-------------------------------
<a name="bond_type"/> bond_type
-------------------------------

**Syntax**:

	bond_type	<type-1> <type-2> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines a bond type and its related attributes. At least one bond type must necessarily be defined before any chemical [bond] is created. A bond type is identified by the types of the two atoms involved.

The parameters _type-1_ and _type-2_ are identifiers of previously defined atom types. Either _type-1_ _type-2_ or _type-2_ _type-1_ will result in the same bond type. Wildcard characters (? and *) can be used to refer to multiple atom types. If a type [prefix] has been previously activated, then the actual identifiers will contain such prefix, followed by the strings in _type-1_ and _type-2_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS bond style associated with the bond type in question.

It is possible to define multiple identical bond types. This can be done, for instance, when the energy and force due to bond stretching is modeled as a sum of multiple terms, each one corresponding to a different model in LAMMPS. Note that bond types can be identical, but cannot be ambiguous.

**Examples**:

	bond_type	CH* CH* 95.877 1.54

In the example above, the given attributes will be associated with every bond involving atoms whose types have names starting with _CH_.

**See also**:

[atom_type], [bond], [write], [prefix]

---------------------------------
<a name="angle_type"/> angle_type
---------------------------------

**Syntax**:

	angle_type	<type-1> <type-2> <type-3> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines an angle type and its related attributes. An angle type is identified by the types of the three atoms involved.

The parameters _type-1_, _type-2_, and _type-3_ are identifiers of previously defined atom types. Either _type-1_ _type-2_ _type-3_ or _type-3_ _type-2_ _type-1_ will result in the same angle type. Wildcard characters (? and *) can be used to refer to multiple atom types. If a type [prefix] has been previously activated, then the actual identifiers will contain such prefix, followed by the strings in _type-1_, _type-2_, and _type-3_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS angle style associated with the angle type in question.

It is possible to define multiple identical angle types. This can be done, for instance, when the energy and force due to a angle bending is modeled as a sum of multiple terms, each one corresponding to a different model in LAMMPS. Note that angle types can be identical, but cannot be ambiguous.

Playmol detects angles automatically as chemical bonds are created. If a detected angle fits to a previously defined angle type, then Playmol will add it to a list of angles to be used later on, for instance, in the creation of a LAMMPS configuration file. If the corresponding angle type does not exist at the moment a new angle is detected, then Playmol will ignore it and produce a warning message.

**Examples**:

	angle_type	CH3 CH2 CH2 62.0965 114

**See also**:

[atom_type], [write], [prefix]

---------------------------------------
<a name="dihedral_type"/> dihedral_type
---------------------------------------

**Syntax**:

	dihedral_type	<type-1> <type-2> <type-3> <type-4> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines a dihedral type and its related attributes. A dihedral type is identified by the types of the four atoms involved.

The parameters _type-1_, _type-2_, _type_3, and _type-4_ are identifiers of previously defined atom types. The order of such identifiers is relevant in the case of dihedral types. Wildcard characters (? and *) can be used to refer to multiple atom types. If a type [prefix] has been previously activated, then the actual identifiers will contain such prefix, followed by the strings in _type-1_, _type-2_, and _type-3_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS dihedral style associated with the dihedral type in question.

It is possible to define multiple identical dihedral types. This can be done, for instance, when the energy and force due to a dihedral angle bending (i.e. bond torsion) is modeled as a sum of multiple terms, each one corresponding to a different model in LAMMPS. Note that dihedral types can be identical, but cannot be ambiguous.

Playmol detects dihedrals automatically as chemical bonds are created. If a detected dihedral fits to a previously defined dihedral type, then Playmol will add it to a list of dihedrals to be used later on, for instance, in the creation of a LAMMPS configuration file. If the corresponding dihedral type does not exist at the moment a new dihedral is detected, then Playmol will ignore it and produce a warning message.

**Examples**:

	dihedral_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0

**See also**:

[atom_type], [write], [prefix]

---------------------------------------
<a name="improper_type"/> improper_type
---------------------------------------

**Syntax**:

	improper_type	<type-1> <type-2> <type-3> <type-4> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines an improper type and its related attributes. An improper type is identified by the types of the four atoms involved.

The parameters _type-1_, _type-2_, _type_3, and _type-4_ are identifiers of previously defined atom types. The order of such identifiers is relevant in the case of improper types. Wildcard characters (? and *) can be used to refer to multiple atom types. If a type [prefix] has been previously activated, then the actual identifiers will contain such prefix, followed by the strings in _type-1_, _type-2_, and _type-3_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS improper style associated with the improper type in question.

It is possible to define multiple identical improper types. This can be done, for instance, when the energy and force due to a improper angle bending (i.e. bond torsion) is modeled as a sum of multiple terms, each one corresponding to a different model in LAMMPS. Note that improper types can be identical, but cannot be ambiguous.

Impropers are created manually using the command [improper], which requires a previously defined improper type.

**Examples**:

	improper_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0

**See also**:

[atom_type], [improper], [write], [prefix]

---------------------
<a name="atom"/> atom
---------------------

**Syntax**:

	atom		<id> <type>

* _id_ = a name to be assigned to the atom being created
* _type_ = the name of a previously defined atom type

**Description**:

This command creates an atom of a given type.

The parameter _id_ must be a single string with no comment tags (#) and no wildcard characters (* or ?). If an atom [prefix] has been previously activated, then the actual atom identifier will contain such prefix, followed by the text specified as _id_. Two atoms cannot have the same identifier.

The parameter _type_ is the identifier of a previously defined atom type. A unique identifier must be provided, with no use of wildcard characters (* or ?). If a type [prefix] has been previously activated, then the actual identifier will contain such prefix, followed by the string in _type_. The specified atom type must have a [mass] associated with it.

When an atom is created, a new monoatomic molecule is automatically detected. If there are, for instance, _n_ already existing molecules, the new monoatomic molecule will receive an index _n+1_. If a [bond] is created afterwards connecting this atom to another one, the two involved molecules will be fused together.

**Examples**:

	atom	C1 C
	atom	C2 C
	atom    H1 H

**See also**:

[atom_type], [mass], [bond]

-------------------------
<a name="charge"/> charge
-------------------------

**Syntax**:

	charge		<atom> <value>

* _atom_ = the name of an atom
* _value_ = charge of the specified atom

**Description**:

This command defines the electric charge of an atom. Wildcard characters (* or ?) can be used to assign the same charge value to different atoms. The command applies for atoms defined either beforehand or afterwards.

If a given [atom] has no charge explicitly defined, then it is considered to be neutral.

**Examples**:

	charge		O* -1.0

In the example above, a change  value of _-1.0_ is assigned to all atoms whose identifiers start with O.

**See also**:

[atom]

---------------------
<a name="bond"/> bond
---------------------

**Syntax**:

	bond		<atom-1> <atom-2>

* _atom-x_ = the name of a previously defined atom

**Description**:

This command creates a chemical bond between two atoms.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be provided, with no use of wildcard characters (* or ?).

If a bond is created between atoms that belong to the same molecule, then a cyclic structure is formed and the number of molecules remains the same. However, if the two atoms belong to distinc molecules, then these molecules are fused together and the number of molecules is decreased. Supposing that the fused molecules have indices _i_ and _j_, the index of the new molecule will be _min(i,j)_. In order to keep a consistent sequence, all existing molecules with indices greater than _max(i,j)_ will have their indices decreased in one unit.

New angles and new dihedrals are detected automatically when a bond is created. If a detected angle/dihedral fits to a previously defined type, then Playmol will add it to a list of angles/dihedrals to be used later on, for instance, in the creation of a LAMMPS configuration file. If the corresponding type does not exist at the moment a new angle/dihedral is detected, then Playmol will ignore it and produce a warning message.

At any moment after chemical bonds have been created, one can check the current list of molecules by invoking the command [write] with the _summary_ format.

**Examples**:

	bond	C1 C2

**See also**:

[bond_type], [atom_type]

-----------------------------
<a name="improper"/> improper
-----------------------------

**Syntax**:

	improper	<atom-1> <atom-2> <atom-3> <atom-4>

* _atom-x_ = name of a previously defined atom

**Description**:

This command creates an improper involving the specified atoms.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be provided, with no use of wildcard characters (* or ?). A previously defined improper type is required, noting that the order of the involved atom types is relevant.

IMPORTANT: Impropers are not detected automatically because there are several possible improper definitions.

**Examples**:

	improper	C1 C2 C3 C4

**See also**:

[improper_type], [atom]

-----------------------------------------
<a name="extra_dihedral"/> extra_dihedral
-----------------------------------------

**Syntax**:

	extra_dihedral	<atom-1> <atom-2> <atom-3> <atom-4>

* _atom-x_ = name of a previously defined atom

**Description**:

This command creates an extra dihedral involving the specified atoms.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be provided, with no use of wildcard characters (* or ?). A previously defined dihedral type is required.

IMPORTANT: Extra dihedrals are only required by some molecular models that define unconventional dihedrals, for example, involving non-bonded atoms. Most models, however, will involve only automatically detectable dihedrals.

**Examples**:

	extra_dihedral	C1 C2 C3 C4

**See also**:

[dihedral_type], [atom]

-------------------
<a name="xyz"/> xyz
-------------------

**Syntax**:

	xyz     	[<file>]

* _file_ (optional) = name of a file containing atom coordinates

**Description**:

This command specifies atomic coordinates for one or more molecules.

The parameter _file_ is the name of a text file containing atomic coordinates. The required format is the [xyz file format], but with element symbols replaced by [atom] identifiers. The parameter _file_ is optional. If it is omitted, then the coordinates must be provided in the subsequent lines, also adhering to the xyz format. In this sense, it becomes completely equivalent to invoke _xyz <file>_ or to invoke _xyz_ followed by _include <file>_.

The xyz format for Playmol is:

* The first line contains the number _N_ of coordinates to be provided
* The second line contains a comment or title
* The _N_ subsequent lines contain four fields separated be spaces or tabs: <atom-id>  <x>  <y>  <z>

The provided coordinates are meant to define a list of actual structures for the molecules previously assembled using the commands [atom] and [bond]. Multiple molecules can be specified either using separate _xyz_ commands or using a unique _xyz_ command involving all coordinates together. The coordinates of all atoms of a given molecule must be provided in sequence, but not in any specific order. It is possible to specify many copies of the same molecule, but a given atom cannot reappear until all other atoms of the same molecule have been defined.

The list of molecular structures created via _xyz_ will be employed to define a simulation box if the command [write] is invoked later on. The origin of the Cartesian space will be located in the center of the simulation box.

**Examples**:

	xyz	H2O.xyz

In the example above, the atomic coordinates of one or more copies of a water molecule is read from the file _H2O.xyz_.

	xyz

	6
	# TIP3P water molecules:
	O   0.00000  0.00000 -1.50000
	H1 -0.58588  0.75695 -1.50000
	H2  0.58588  0.75695 -1.50000
	O   0.00000  0.00000  1.50000
	H1 -0.58588 -0.75695  1.50000
	H2  0.58588 -0.75695  1.50000

In the example above, the atomic coordinates of two water molecules are provided inline in the Playmol script.

**See also**:

[atom], [bond], [box], [write], [include]

-------------------
<a name="box"/> box
-------------------

**Syntax**:

	box		lengths <lx> <ly> <lz>
	box		volume <value> [aspect <ax> <ay> <az>]
	box		density <value> [aspect <ax> <ay> <az>]

* _lx_, _ly_, _lz_ = box dimensions in directions x, y, and z
* _value_ = density of the simulation box
* _ax_, _ay_, _az_ (optional) = Scaling factors in directions x, y, and z

**Description**:

This command specifies the side lengths, the volume, or the density of a simulation box.

The parameters _lx_, _ly_, and _lz_ are the dimensions of the simulation box in the three spatial directions.

The parameter _value_ is either the volume or the density of the simulation box, depending on the kind of specification being done.

The optional keyword _aspect_ with parameters _ax_, _ay_, and _az_ can be used to determine the relative box sizes in the three spatial directions (the absolute values of these three parameters are meaningless). If the keyword _aspect_ is omitted, Playmol will consider a cubic box.

IMPORTANT: when the box density is specified, Playmol will automatically calculate the box side lengths when necessary. For this, it will sweep the whole list of molecular structures defined via [xyz] commands in order to determine the total mass of the system. Since Playmol will not perform unit conversion, the specified density value must be consistent with the mass and length units used in other commands. For instance, if [LAMMPS real units] are considered (mass in g/mol and lengths in Å), then the density must be provided in Da/Å³ (therefore, values originally in g/cm³ must be multiplied by 0.60221413).

When building a simulation box, Playmol places its geometric center in the origin of the Cartesian space. To allow periodic boundary conditions, Playmol does not check whether the previously specified atomic coordinates fit inside the central box.

**Examples**:

	box	size 30 30 30

The example above specifies a cubic simulation box.

	box	volume 27000

The example above is completely equivalent to the preceding one, since the default aspect is cubic.

	box	density 0.602214 aspect 1 1 2

The example above specifies a simulation box with density equal to 0.602214 in LAMMPS "real" units (equivalent to 1.0 g/cm³). The box is elongated twice in the z direction with respect to the other two directions.

**See also**:

[xyz], [packmol], [write]

---------------------------
<a name="packmol"/> packmol
---------------------------

**Syntax**:

	packmol		<keyword> <arguments> [<keyword> <arguments> ...]

* _keyword_ =  _seed_ or _tolerance_ or _nloops_ or _retry_ or _fix_ or _copy_ or _pack_ or _action_

**Keywords and arguments**:

* seed <iseed>
	* _iseed_ = an integer number for seeding Packmol's random number generator (_default_ = 1234)
* tolerance <tol>
	* _tol_ = the minimum distance to be observed between atoms from distinct molecules (_default_ = 1.0)
* nloops <N>
	* _N_ = the maximum number of iterations of the Packmol algorithm (_default_ = 50)
* retry <factor>
	* _factor_ = a scaling factor for tolerance values in successive packing attempts (_default_ = 1.0)
* fix <molecule> <x> <y> <z>
	* _molecule_ = the index of an existing molecule with predefined coordinates
	* _x_, _y_, _z_ = a spatial coordinate for placing the molecule's geometric center
* copy <molecule> <N>
	* _molecule_ = the index of an existing molecule with predefined coordinates
	* _N_ = the number of randomly positioned copies of the molecule
* pack <molecule> <N>
	* _molecule_ = the index of an existing molecule with predefined coordinates
	* _N_ = the number of randomly packed copies of the molecule
* action <option>
	* _option_ = _execute_ or _setup_

**Description**:

This command creates Packmol input files or invokes Packmol to build a low-overlap molecular packing inside a previously specified simulation box.

The keywords _seed_, _tolerance_, and _nloops_ change the values of some parameters that affect Packmol's behavior. Different integer values for _seed_ tell Packmol to use different random number sequences for its packing algorithm. The parameter _tolerance_ is the minimum distance that atoms of distinct molecules must keep from each other in the final packing. The parameter _nloops_ is the maximum number of iterations to be carried out with the Packmol algorithm in a packing attempt.

The parameter _retry_ is a reduction factor. If its value is 1.0, Playmol will invoke Packmol only once with the specified tolerance and will produce a warning message if the packing fails. If _retry_ is smaller than 1.0, then Playmol will invoke Packmol as many times as necessary to achieve a successful packing, with tolerance being iteratively multiplied by the value of _retry_ at each new attempt.

The keywords _fix_, _copy_, and _pack_ require the index of an existing molecule for which at least one set of atomic coordinates have been defined using the command [xyz]. The first defined set of coordinates for such molecule will then be used as a rigid-body model for translation or replication. The meaning of each option is:

* fix <index> <x> <y> <z>: makes a copy of Molecule _index_ with its geometric center located at the provided coordinate (_x_, _y_, _z_), keeping its original orientation.

* copy <index> <N>: makes _N_ copies of Molecule _index_ in random positions, but keeping the original orientation (i.e. all _N_ copies will be aligned in the final packing). This is useful for packing long molecules.

* pack <index> <N>: makes _N_ copies of Molecule _index_ in random positions and with random orientations.

IMPORTANT: the option _summary_ of command [write] can be very useful for checking out the indices of the molecules.

The keyword _action_ is used to create Packmol input files or to invoke Packmol. The following options are available:

* __execute__: this option calls Packmol to build the desired molecular packing. It requires the previous definition of a simulation [box]. Moreover, it requires that at least one keyword _fix_, _copy_, or _pack_ has appeared in a previous packmol command or appears in the same packmol command, either before or after the keyword _action_. If the parameter _retry_ is currently equal to 1.0 (its default value), then Packmol will do only one packing attempt with the specified _tolerance_ and produce a warning message in case such attempt fails. If _retry_ is smaller than 1.0, then Packmol will keep trying until a successful attempt is achieved, with _tolerance_ being iteratively multiplied by the _retry_ value at each new attempt. IMPORTANT: if packing succeeds, then the current list of atomic coordinates is replaced by the new coordinates generated by Packmol. After that, the command [write] can be used to create a LAMMPS configuration file with the attained packing.

* __setup__: this option generates an input file named _packmol.inp_ and coordinate files _molecule-x.inp_, where _x_ is the index of each involved molecule. These files are prepared for running Packmol externally in order to generate a file _packmol-output.xyz_ containing the final packing if the algorithm succeeds with the given tolerance. For illustration, one may notice that a successful use of the packmol command with options _retry 1.0 action execute_ would have exactly the same result as the following sequence of commands:


	packmol 	action setup
	shell   	packmol < packmol.inp
	reset		xyz
	xyz     	packmol_output.xyz

Nevertheless, the real usefulness of the option _setup_ is to permit editing of the file _packmol.inp_ in order to impose some additional constraints, which are not directly handled in the current version of Playmol. Please see [Packmol User's Guide] for additional information.

The Packmol algorithm is described in the following paper:

	@article{packmol_paper,
		author = {Martinez, L. and Andrade, R. and Birgin, E. G. and Martinez, J. M.},
		title = {PACKMOL: A package for building initial configurations for molecular dynamics simulations},
		journal = {J. Comput. Chem.},
		volume = {30},
		number = {13},
		pages = {2157--2164},
		year = {2009}
		doi = {10.1002/jcc.21224},
	}

**Examples**:

	box    		density 0.602214
	packmol		tolerance 3.0 retry 0.9 fix 1 0.0 0.0 0.0 pack 2 1000 action execute
	write		lammps system.data

The example above uses Packmol to create a random packing of molecules with density equal to 0.602214 Da/Å³ (1.0 g/cm³) in which one copy of molecule 1 is centered at the origin and 1000 copies of molecule 2 are packed with random positions and random orientations. The desired mininum intermolecular atomic distance is initially set to 3.0 Å, but Packmol will keep reducing the tolerance in 90% until the packing is successful. Finally, a LAMMPS configuration file is generated.

**See also**:

[box], [xyz], [write]

-----------------------
<a name="write"/> write
-----------------------

**Syntax**:

	write		 <format> [<file>]

* _format_ = _playmol_ or _lammps_ or _summary_ or _xyz_ or _lammpstrj_
* _file_ (optional) = name of a file to be created

**Description**:

This command writes down the molecular system in a file or in the standard output device.

The parameter _format_ must be one of the following options:

* __playmol__: the output will contain Playmol commands that could be used in another script to build the same system. For illustration, detected angles and dihedrals appear as commented lines. Type and atom prefixes are explicitly added to the corresponding identifiers.

* __lammps__: the command will produce information in the LAMMPS configuration file format, which can be used as an initial configuration for a Molecular Dynamics simulation using LAMMPS' command [read_data].

* __summary__: this option will print a summary of the system characteristics, including the amount of every defined and detected structure such as angles, dihedrals, and molecules. Properties of each molecular species will also be printed, such as its mass, its atoms, and the number of defined sets of coordinates. This is useful for debugging purposes.

* __xyz__: writes down the list of atomic coordinates using the [xyz file format], but with element symbols replaced by atom identifiers. This is useful for using with another Playmol script or for visualization purposes.

* __lammpstrj__: writes down the list of atomic coordinates using the LAMMPS trajectory format. This is useful for visualization with [VMD].

The optional parameter _file_ is the name of the file which will contain the system info. If it is omitted, the info will be written in the standard output unit (the computer screen, in most cases).

**Examples**:

	write	playmol water.mol
	write	lammps water.data

-------------------------
<a name="prefix"/> prefix
-------------------------

**Syntax**:

	prefix		<target> <string>

* _target_ = _types_ or _atoms_
* _string_ = _none_ or a character string to be used as prefix for atom types or atoms

**Description**:

This command enables or disables a prefix to be added to every atom type identifier or to every atom identifier in subsequent commands.

The parameter _target_ indicates which type of identifier will be modified by the enabled prefix. Use _types_ to apply the prefix to atom type identifiers or _atoms_ to apply the prefix to atom identifiers.

The parameter _string_ must be a single character string with no comment tags (#) and no wildcard characters (* or ?). Such string will be used in subsequent commands as a prefix for atom type identifiers or for atom identifiers, depending on the specified parameter _target_. If _string_ is _none_, then no prefix will be used for the specified _target_ in subsequent commands.

Both prefix and [suffix] can be defined simultaneously.

**Examples**:

	prefix		atoms H2O-
	atom		H1 H
	atom		H2 H

The example above creates two atoms whose identifiers are H2O-H1 and H2O-H2.

**See also**:

[atom_type], [atom], [suffix]

-------------------------
<a name="suffix"/> suffix
-------------------------

**Syntax**:

	suffix		<target> <string>

* _target_ = _types_ or _atoms_
* _string_ = _none_ or a character string to be used as suffix for atom types or atoms

**Description**:

This command enables or disables a suffix to be added to every atom type identifier or to every atom identifier in subsequent commands.

The parameter _target_ indicates which type of identifier will be modified by the enabled suffix. Use _types_ to apply the suffix to atom type identifiers or _atoms_ to apply the suffix to atom identifiers.

The parameter _string_ must be a single character string with no comment tags (#) and no wildcard characters (* or ?). Such string will be used in subsequent commands as a suffix for atom type identifiers or for atom identifiers, depending on the specified parameter _target_. If _string_ is _none_, then no suffix will be used for the specified _target_ in subsequent commands.

Both [prefix] and suffix can be defined simultaneously.

**Examples**:

	suffix		atoms -H2O
	atom		H1 H
	atom		H2 H

The example above creates two atoms whose identifiers are H1-H2O and H2-H2O.

**See also**:

[atom_type], [atom], [prefix]

---------------------------
<a name="include"/> include
---------------------------

**Syntax**:

	include		<file>

* _file_ = name of a file containing additional playmol commands

**Description**:

This command redirects Playmol to read and execute commands from another script.

**Examples**:

	include		H2O.mol

**See also**:

[atom_type], [atom]

-----------------------
<a name="reset"/> reset
-----------------------

**Syntax**:

	reset		<list>

* _list_ = *bond_types* or *angle_types* or *dihedral_types* or *improper_types* or *atoms* or *charges* or *bonds* or *impropers* or *xyz* or *packmol* or *all*

**Description**:

This command resets one or more lists of predefined entities.

The options *bond_types*, *angle_types*, *dihedral_types*, *improper_types*, *charges*, and *impropers* are self-explanatory. The remaining options are:

* _atoms_: resets the list of atoms and its dependent lists: charges, bonds, angles, dihedrals, impropers, molecules, coordinates, and Packmol definitions.

* _bonds_: resets the list of bonds and all its dependent lists: angles, dihedrals, molecules, coordinates, and Packmol definitions.

* _xyz_: resets the list of coordinates.

* _packmol_: resets the list of Packmol definitions.

* _all_: resets all lists, including atom types and masses.

**See also**:

[packmol]

-----------------------
<a name="shell"/> shell
-----------------------

**Syntax**:

	shell		<command>

* _command_ = an external shell command

**Description**:

This command executes an external shell command.

**Example**:

	shell	mv packmol.inp system.inp
	shell	packmol < system.inp

The example above uses Linux command _mv_ to rename a file from _packmol.inp_ to _system.inp_ and then executes Packmol using it as input.

**See also**:

[packmol]

---------------------
<a name="quit"/> quit
---------------------

**Syntax**:

	quit [all]

**Description**:

This command interrupts the execution of a Playmol script, thus being useful for debugging purposes.

The optional keyword _all_ forces Playmol to stop completely. In the absence of the keyword, the command _quit_ has the same behavior of a normal end of file. For instance, if the script has been invoked by another script via the [include] command, execution of the invoking script continues.

**Example**:

	write summary
	quit

The example above writes a summary of the current molecular system and then quits Playmol.

**See also**:

[include]

---------------------

[atom_type]:		#atom_type
[mass]:			#mass
[bond_type]:		#bond_type
[angle_type]:		#angle_type
[dihedral_type]:	#dihedral_type
[improper_type]:	#improper_type
[atom]:			#atom
[charge]:		#charge
[bond]:			#bond
[extra_dihedral]:	#extra_dihedral
[improper]:		#improper
[xyz]:			#xyz
[box]:			#box
[packmol]:		#packmol
[write]:		#write
[prefix]:		#prefix
[suffix]:		#suffix
[include]:		#include
[reset]:		#reset
[shell]:		#shell
[quit]:			#quit

[LAMMPS real units]:	http://lammps.sandia.gov/doc/units.html
[read_data]:		http://lammps.sandia.gov/doc/read_data.html
[xyz file format]:	http://openbabel.org/wiki/XYZ_(format)
[Packmol User's Guide]:	http://www.ime.unicamp.br/~martinez/packmol/quickguide/
[VMD]:			http://www.ks.uiuc.edu/Research/vmd/

