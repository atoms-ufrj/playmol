Commands      {#commands}
========

Playmol is designed to execute text-file scripts containing the commands described below:

* [atom_type] - creates an atom type with given name and parameters.
* [mass] - specifies the mass of atoms of a given type.
* [bond_type] - defines parameters of bonds between atoms of two given types.
* [angle_type] - defines parameters of angles involving atoms of three given types.
* [dihedral_type] - defines parameters of dihedrals involving atoms of four given types.
* [improper_type] - defines parameters of impropers involving atoms of four given types.
* [atom] - creates an atoms with given name and type.
* [charge] - specifies the charge of a given atom.
* [bond] - creates a bond between two given atoms (angles and dihedrals are automatically detected).
* [improper] - creates an improper involving a given atom quadruplet.
* [xyz] - defines positions for all atoms of one or more molecules.
* [box] - defines the properties of a simulation box.
* [packmol] - generates Packmol input files or execute Packmol to create a packed molecular system.
* [write] - writes down system info in different file formats (including _LAMMPS_ data files).
* [prefix] - defines default prefixes for atom types and atoms.
* [include] - includes all commands of another script.
* [quit] - interrupts Playmol execution.

[atom_type]:		#atom_type
[mass]:			#mass
[bond_type]:		#bond_type
[angle_type]:		#angle_type
[dihedral_type]:	#dihedral_type
[improper_type]:	#improper_type
[atom]:			#atom
[charge]:		#charge
[bond]:			#bond
[improper]:		#improper
[xyz]:			#zyx
[box]:			#box
[packmol]:		#packmol
[write]:		#write
[prefix]:		#prefix
[include]:		##include
[quit]:			#quit

-------------------------------
<a name="atom_type"/> atom_type
-------------------------------

**Syntax**:

    atom_type	id attribute-list

* _id_ = a name assigned to the atom type being defined
* _attribute-list_ = a list of attributes related to the atom type being defined

**Examples**:

    atom_type	CH3 0.1947 3.75
    atom_type	CH2 lj/cut 0.0914 3.95

**Description**:

This command defines an atom type and its related attributes. At least one atom type must necessarily be defined before any [atom] is created and before any [bond_type], [angle_type], [dihedral_type], or [improper_type] is
defined.

The parameter _id_ must be a single string with no comment tags (#) and no wildcard characters (* or ?). If a type [prefix] has been previously activated, then the actual atom type identifier will contain such prefix, followed by the text specified as _id_. Two atom types cannot have the same identifier.

The parameter _attribute-list_ can be a sequence of strings containing any characters except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS pair style associated with the atom type in question.

**See also**:

[atom], [bond_type], [angle_type], [dihedral_type], [improper_type], [prefix]

---------------------
<a name="mass"/> mass
---------------------

**Syntax**:

    mass		type value

* _type_ = name of a previously defined atom type
* _value_ = mass of atoms of the specified atom type

**Examples**:

    mass		CH2 14.0
    mass		H* 1.0

**Description**:

This command defines the mass of all atoms of a given type.

**See also**:

[atom_type]

-------------------------------
<a name="bond_type"/> bond_type
-------------------------------

**Syntax**:

    bond_type	type-1 type-2 attribute-list

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Examples**:

    bond_type	CH3 CH2 95.877 1.54
    bond_type	CH2 CH? 95.877 1.54
    bond_type	C* C* 95.877 1.54

**Description**:

This command defines a bond type and its related attributes. One or more bond types must necessarily be defined before any chemical [bond] is created. A bond type is identified by the types of the two atoms involved.

The parameters _type-1_ and _type-2_ are identifiers of previously defined atom types. Wildcard characters (? and *) can be used to refer to multiple atom types. Two bond types cannot refer to the same pair of atom types, regardless of the order of _type-1_ and _type-2_. If a type [prefix] has been previously activated, then the actual identifiers of both atom types will contain such prefix, followed by the specified parameters _type-1_ and _type-2_.

The parameter _attribute-list_ can be a sequence of strings containing any characters except comment tags (#). When a LAMMPS data file is generated using the command [write], then _attribute-list_ is assumed to contain the parameters of the LAMMPS bond style associated with the bond type in question.

**See also**:

[atom_type], [bond], [write]

---------------------------------
<a name="angle_type"/> angle_type
---------------------------------

**Syntax**:

    angle_type	atom-type-1 atom-type-2 atom-type-3 attribute-list

* _atom-type-x_ = name of a previously defined atom type
* _attribute-list_ = list of attributes related to the defined angle type

**Examples**:

    angle_type	CH3 CH2 CH2 62.0965 114
    angle_type	CH2 CH2 CH2 62.0965 114

**Description**:

This command defines....

**See also**:

[atom_type]

---------------------------------------
<a name="dihedral_type"/> dihedral_type
---------------------------------------

**Syntax**:

    dihedral_type	atom-type-1 atom-type-2 atom-type-3 atom-type-4 attribute-list

* atom-type-x = name of a previously defined atom type
* attribute-list = list of attributes related to the defined dihedral type

**Examples**:

    dihedral_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0
    dihedral_type	CH2 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0

**Description**:

This command defines....

**See also**:

[atom_type]

---------------------------------------
<a name="improper_type"/> improper_type
---------------------------------------

**Syntax**:

    improper_type	atom-type atom-type atom-type atom-type attribute-list

* _atom-type-x_ = name of a previously defined atom type
* _attribute-list_ = list of attributes related to the defined improper type

**Examples**:

    improper_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0
    improper_type	CH2 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0

**Description**:

This command defines....

**See also**:

[atom_type], [improper]

---------------------
<a name="atom"/> atom
---------------------

**Syntax**:

    atom		name atom-type

* _name_ = name of the atom being defined
* _atom-type_ = name of a previously defined atom type

**Examples**:

    atom		C1 CH3
    atom		C2 CH2

**Description**:

This command defines.....

Both parameters @a name and @a atom-type must not contain wildcard characters.

Everytime a new atom is created, a new molecule is also defined which contains such atom. Afterwards,
if the new atom is bonded to another existing atom via the [bond] command, the two corresponding
molecules will be fused together.

**See also**:

[atom_type], [bond]

---------------------
<a name="bond"/> bond
---------------------

**Syntax**:

    bond		atom-1 atom-2

* _atom-x_ = name of a previously defined atom

**Examples**:

    bond		C1 C2

**Description**:

This command defines....

**See also**:

[bond_type], [atom_type]

-----------------------------
<a name="improper"/> improper
-----------------------------

**Syntax**:

    improper	atom-1 atom-2 atom-3 atom-4

* _atom-x_ = name of a previously defined atom

**Examples**:

    improper	C1 C2 C3 C4

**Description**:

This command defines....

**See also**:

[improper_type], [atom]

-------------------
<a name="xyz"/> xyz
-------------------

**Syntax**:

    xyz		[file-name]

* _file-name_ (optional) = name of a file containing atom coordinates

**Examples**:

    xyz		H2O.xyz

    xyz
    3
    # Water molecule:
    H1		-0.239 0.928 0.000
    O		 0.000 0.000 0.000
    H2		 0.958 0.000 0.000

**Description**:

This command.....

**See also**:

[atom], [box]

-------------------
<a name="box"/> box
-------------------

**Syntax**:

    box		x-min x-max y-min y-max z-min z-max

* _x-min_ _x-max_ = simulation box limits in the x-direction
* _y-min_ _y-max_ = simulation box limits in the y-direction
* _z-min_ _z-max_ = simulation box limits in the z-direction

**Examples**:

    box		0.0 30.0 0.0 30.0 0.0 30.0

**Description**:

This command.....

**See also**:
[xyz]

---------------------------
<a name="packmol"/> packmol
---------------------------

**Syntax**:

    packmol tolerance seed molecule id mode [molecule id mode] action

* _tolerance_ =
* _seed_ = 

**Examples**:

    packmol		

**Description**:

This command...

**See also**:

[box], [xyz]

-----------------------
<a name="write"/> write
-----------------------

**Syntax**:

    write		style [file-name]

* _style_ = _playmol_ or _lammps_ or _summary_
* _file-name_ (optional) = name of a file to be created

**Examples**:

    write		playmol water.pmol
    write		lammps water.lammps_data

**Description**:

This command writes down the molecular system in a file (if specified) or in the standard output device.

-------------------------
<a name="prefix"/> prefix
-------------------------

**Syntax**:

    prefix		prefix-style text


* _prefix-style_ = @a type or @a atom
* _text_ = 

**Examples**:

    prefix		type A
    prefix		atom H2O_

**Description**:

This command...

**See also**:

[atom_type], [atom]

---------------------------
<a name="include"/> include
---------------------------

**Syntax**:

    include		file-name

* _file-name_ = name of a file containing other playmol commands

**Examples**:

    include		H2O.pmol

**Description**:

This command.....

**See also**:

[atom_type], [atom]

---------------------
<a name="quit"/> quit
---------------------

**Syntax**:

    quit

**Description**:

This command.....

