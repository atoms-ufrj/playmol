Playmol Commands      {#commands}
================

Playmol is designed to execute scripts containing the commands described below. Click in the command
name for detailed description and examples. In all examples, the units of measurement employed for
physically meaningful values are those corresponding to [LAMMPS real units].

| Command         | Description                                                               |
|:---------------:|:--------------------------------------------------------------------------|
| [define]        | defines a variable for further substitution                               |
| [for/next]      | executes commands repeatedly while changing the value of a variable       |
| [if/then/else]  | executes commands conditionally or selects between two command sequences  |
| [atom_type]     | creates an atom type with given name and parameters                       |
| [mass]          | specifies the mass of atoms of a given type                               |
| [diameter]      | specifies the diameter of atoms of a given type                           |
| [bond_type]     | defines parameters for bonds between atoms of two given types             |
| [angle_type]    | defines parameters for angles involving atoms of three given types        |
| [dihedral_type] | defines parameters for dihedrals involving atoms of four given types      |
| [improper_type] | defines parameters for impropers involving atoms of four given types      |
| [atom]          | creates an atom with given name and type                                  |
| [charge]        | specifies the charge of a given atom                                      |
| [bond]          | creates chemical bonds and automatically detects angles and dihedrals     |
| [improper]      | creates an improper involving four given atoms or search for impropers    |
| [extra]         | creates an extra bond, angle, or dihedral involving given atoms           |
| [link]          | links two atoms (and fuses their molecules) without actually bonding them |
| [unlink]        | removes an existing link (and splits the corresponding molecule)          |
| [build]         | guesses atom positions from provided geometric information                |
| [prefix/suffix] | defines default prefixes or suffixes for atom types and atoms             |
| [box]           | defines the properties of a simulation box                                |
| [align]         | aligns the principal axes of a molecule to the Cartesian axes             |
| [packmol]       | executes Packmol in order to create a packed molecular system             |
| [write]         | saves system info in different file formats (including LAMMPS data files) |
| [include]       | includes commands from another script                                     |
| [reset]         | resets a list of entities together with its dependent lists               |
| [shell]         | executes an external shell command                                        |
| [quit]          | interrupts the execution of a Playmol script                              |

----------------------------------------------------------------------------------------------------
<a name="define"></a>
define
----------------------------------------------------------------------------------------------------

**Syntax**:

	define 	<variable> as <string>

* _variable_ = a valid name to be assigned to the defined variable
* _string_ = a single character string

**Description**:

This command defines a variable whose value (a character string) is intended to be substituted in
forthcoming commands.

The parameter _variable_ the must be a single character string containing only letters (A-Z or a-z),
numbers (0-9), and underscores (_). The first character must necessarily be a letter. __Important__:
variable names are case-sensitive.

The parameter _value_ must be a single character string with no comment tags (#).

Substitution is instructed in forthcoming commands by either preceding the variable name with a
symbol "$" or enclosing it between symbols "${" and "}", with no spacing. Every time the symbol "$"
is found in a command, Playmol will admit that a variable is being referred to. The actual command
will only be issued after all variables have been replaced by their values.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define		name as C
atom_type	$name	0.0914 3.95
atom		${name}1 $name
atom		${name}2 $name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, an atom type named _C_ is defined with parameters 0.0914 and 3.95. Then, two
atoms named _C1_ and _C2_ are created, both of type _C_.

**See also**:

[for/next], [if/then/else], [atom_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="for_next"></a>
for/next
----------------------------------------------------------------------------------------------------

**Syntax**:

	for <variable> from <min> to <max>
		commands
	next

or

	for <variable> from <max> downto <min>
		commands
	next

or

	for <variable> in <string-list>
		commands
	next

* _variable_ = a valid name to be assigned to the defined variable
* _min_, _max_ = either two integer numbers or two single characters, with _max_ &ge; _min_
* _string-list_ = a list of character strings

**Description**:

This command repeatedly executes a sequence of command-lines between the _for_ and _next_ statements
while the value of a variable changes in accordance with a specified range or list.

The parameter _variable_ the must be a valid variable name, that is, the first character must be a
letter (A-Z or a-z) and the remaining characters can be letters (A-Z or a-z), numbers (0-9) or
underscores (_).

The parameters _min_ and _max_ can be either two integer numbers with _max_ &ge; _min_ or two single
characters with _max_ not preceding _min_ in the [ASCII] character table. If the corresponding
condition is not met, then the command-lines between _for_ and _next_ will be ignored. Otherwise,
these command-lines will be executed repeatedly, with the value of _variable_ ranging either from
_min_ up to _max_ or from _max_ down to _min_, depending on the chosen syntax (see above).

The parameter _string-list_ is a list of character strings separated by spaces and/or tabs. If the
list is empty, then the command-lines between _for_ and _next_ will be ignored. Otherwise, these
command-lines with be executed repeatedly, with the value of _variable_ turning into the entries of
_string-list_, one at a time and in the order they appear in the list.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atom	C1 CH3
for i from 2 to {$N-1}
	atom	C$i CH2
next
atom	C$N CH3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the units of a united-atom n-alkane molecule are defined
and named as C1, C2, C3, etc. It is assumed that variable N (the number of
carbon atoms) has already been defined.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for mol in water ethanol glycerol
	build $mol.xyz
next
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the files _water.xyz_, _ethanol.xyz_, and _glycerol.xyz_ are employed as
inputs for the [build] command.

**See also**:

[define], [if/then/else]

----------------------------------------------------------------------------------------------------
<a name="if_then_else"></a>
if/then/else
----------------------------------------------------------------------------------------------------

**Syntax**:

	if <condition> then
		commands
	endif

or

	if <condition> then
		commands
	else
		commands
	endif

* _condition_ = 1 or 0 (e.g. the result of a logical expression)

**Description**:

This command conditionally executes a sequence of command-lines or selects between two given
sequences of command-lines.

The parameter _condition_ can be either _1_ (true) or _0_ (false). If any value other than these two
is passed, an error message is produced. In practice, _condition_ might be the result of a variable
substitution or a logical expression evaluation (see the [Playmol Basics] section for details).

If the first syntax above is used (that is, without the _else_ statement), then the command-lines
between _if_ and _endif_ will be executed if _condition_ = _1_ or ignored if _condition_ == 0.
Otherwise, if the second syntax is chosen, then the executed command-lines will either be those
between _if_ and _else_ if _condition_ = 1 or those between _else_ and _endif_ if _condition_ = 0.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nerd force field for normal alkanes:
define		factor as 1.987E-3
atom_type	CH2 {45.8*$factor}  3.930
if {$N == 2} then
  atom_type	CH3 {100.6*$factor} 3.825
else
  if {$N == 3} then
    atom_type	CH3 {102.6*$factor} 3.857
  else
    atom_type	CH3 {104.0*$factor} 3.910
  endif
endif
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the attribute list of the atom type _CH3_ will depend on the number of carbon
atoms (_N_) of the n-alkane molecule. Ethane and propane have their own parameters, while n-butane
and superior alkanes share the same parameters.

**See also**:

[define], [for/next]

----------------------------------------------------------------------------------------------------
<a name="atom_type"></a>
atom_type
----------------------------------------------------------------------------------------------------

**Syntax**:

	atom_type 	<id> <attribute-list>

* _id_ = a name to be assigned to the atom type being defined
* _attribute-list_ = a list of attributes related to the atom type being defined

**Description**:

This command defines an atom type with its related attributes. At least one atom type must
necessarily be defined before any [atom] is created and before any [bond_type], [angle_type],
[dihedral_type], or [improper_type] is defined.

The parameter _id_ must be a single string with no comment tags (#) and no wildcard characters
(* or ?). If a type-related [prefix/suffix] has been previously activated, then the actual atom type
identifier will contain such prefix and/or suffix added to _id_. Distinct atom types cannot have the
same identifier.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except
comment tags (#). When a [LAMMPS data file] is generated using the command [write], then the content
of _attribute-list_ is assumed to be the parameters of the [LAMMPS pair style] associated with the
atom type in question.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atom_type	O	0.1553 3.166
atom_type	H	0.0000 0.000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, atom types _O_ and _H_ are defined with Lennard-Jones parameters _epsilon_ and
_sigma_, respectively, corresponding to the SPC water model. This can be used to generate a LAMMPS
configuration file associated with the pair style _lj/cut/coul/long_, for instance.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prefix		types TraPPE-
atom_type	CH3 lj/cut 0.1947 3.75
atom_type	CH2 lj/cut 0.0914 3.95
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, united-atom types _TraPPE-CH3_ and _TraPPE-CH2_ are defined Lennard-Jones
parameters _epsilon_ and _sigma_, respectively, corresponding to the TraPPE force field. An explicit
pair style specification _lj/cut_ is done, supposing that a pair style _hybrid_ or _hybrid/overlay_
will be employed in LAMMPS.

**See also**:

[atom], [bond_type], [angle_type], [dihedral_type], [improper_type], [prefix/suffix]

----------------------------------------------------------------------------------------------------
<a name="mass"></a>
mass
----------------------------------------------------------------------------------------------------

**Syntax**:

	mass		<type> <value>

* _type_ = the name of a atom type
* _value_ = mass of any atom with the specified atom type

**Description**:

This command defines the mass of atoms of a given type.

The parameter _type_ is an [atom_type] identifier. Wildcard characters (* or ?) can be used so that
the same mass value is assigned to multiple atom types. If a type-related [prefix/suffix] has been
previously activated, then the actual atom type identifier will contain such prefix and/or suffix
added to _type_. Distinct atom types cannot have the same identifier. __Important__: the presence of
wildcards makes a mass value applicable to atom types defined either beforehand or afterwards.

An [atom] will only be created if a mass value has been previously defined to its corresponding atom
type. For models which contain massless atoms, such as a 4-point or a 5-point water model, one must
define their masses as 0 (zero).

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mass		CH2	14.0
mass		H*	1.008
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, a mass value of _14.0_ is assigned to atom type _CH2_ and a mass value of
_1.008_ is assigned to all atoms types whose identifiers start with _H_.

**See also**:

[atom_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="diameter"></a>
diameter
----------------------------------------------------------------------------------------------------

**Syntax**:

	diameter	<type> <value>

* _type_ = the name of a atom type
* _value_ = diameter of any atom with the specified atom type

**Description**:

This command defines the diameter of atoms of a given type.

The parameter _type_ is an [atom_type] identifier. Wildcard characters (* or ?) can be used so that
the same diameter value is assigned to multiple atom types. If a type-related [prefix/suffix] has
been previously activated, then the actual atom type identifier will contain such prefix and/or
suffix added to _type_. Distinct atom types cannot have the same identifier. __Important__: the
presence of wildcards makes a diameter value applicable to atom types defined either beforehand or
afterwards.

It is not necessary to define a diameter for every atom type. This is only required if one wishes to
setup atom-specific distance tolerances when creating a molecular packing with the [packmol]
command, thus superseding its parameter _tolerance_.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diameter	CH2	3.0
diameter	H*	2.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, a diameter value of _3.0_ is assigned to atom type _CH2_ and a diameter value
of _2.0_ is assigned to all atoms types whose identifiers start with _H_.

**See also**:

[atom_type], [atom], [packmol]

----------------------------------------------------------------------------------------------------
<a name="bond_type"></a>
bond_type
----------------------------------------------------------------------------------------------------

**Syntax**:

	bond_type	<type-1> <type-2> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines a bond type and its related attributes. At least one bond type must necessarily
be defined before any chemical [bond] is created. A bond type is identified by the types of the two
atoms involved.

The parameters _type-1_ and _type-2_ are identifiers of previously defined atom types. Either order
_type-1_ _type-2_ or _type-2_ _type-1_ will result in the same bond type. Wildcard characters (? and
*) can be used to refer to multiple atom types. If a type-related [prefix/suffix] has been
previously activated, then each actual atom type identifier will contain such prefix and/or suffix
added to _type-x_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except
comment tags (#). When a [LAMMPS data file] is generated using the command [write], then the content
of _attribute-list_ is assumed to be the parameters of the [LAMMPS bond style] associated with the
bond type in question.

It is possible to define multiple identical bond types. This can be done, for instance, when the
bond stretching potential is modeled as a sum of multiple terms, each one corresponding to a
different [LAMMPS bond style]. Note that two or more bond types can be defined identically, even
with the use of wildcard characters. However, they cannot be defined ambiguously.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bond_type	CH* CH* 95.877 1.54
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the given attributes will be associated with every bond involving atoms whose
types have identifiers starting with _CH_.

**See also**:

[atom_type], [bond], [write], [prefix/suffix]

----------------------------------------------------------------------------------------------------
<a name="angle_type"></a>
angle_type
----------------------------------------------------------------------------------------------------

**Syntax**:

	angle_type	<type-1> <type-2> <type-3> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines an angle type and its related attributes. An angle type is identified by the
types of the three atoms involved.

The parameters _type-1_, _type-2_, and _type-3_ are identifiers of previously defined atom types.
Either order _type-1_ _type-2_ _type-3_ or _type-3_ _type-2_ _type-1_ will result in the same angle
type. Wildcard characters (? and *) can be used to refer to multiple atom types. If a type-related
[prefix/suffix] has been previously activated, then each actual atom type identifier will contain
such prefix and/or suffix added to _type-x_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except
comment tags (#). When a [LAMMPS data file] is generated using the command [write], then the content
of _attribute-list_ is assumed to be the parameters of the [LAMMPS angle style] associated with the
angle type in question.

It is possible to define multiple identical angle types. This can be done, for instance, when the
angle bending potential is modeled as a sum of multiple terms, each one corresponding to a different
[LAMMPS angle style]. Note that two or more angle types can be defined identically, even with the
use of wildcard characters. However, they cannot be defined ambiguously.

Playmol detects angles automatically when chemical bonds are created. If a detected angle fits to a
previously defined angle type, then Playmol will add it to a list of angles to be used later on, for
instance, in the creation of a [LAMMPS data file]. If the corresponding angle type does not exist at
the moment a new angle is detected, then Playmol will ignore such angle and produce a warning
message. Besides the automatic detection, an angle can also be explicitly created using the [extra]
command.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
angle_type	CH3 CH2 CH2 62.0965 114
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the given attributes will be associated with every detected or specified angle
that involves three consecutive atoms whose type identifiers are _CH3_ _CH2_ _CH2_ or _CH2_ _CH2_
_CH3_.

**See also**:

[atom_type], [write], [prefix/suffix], [extra]

----------------------------------------------------------------------------------------------------
<a name="dihedral_type"></a>
dihedral_type
----------------------------------------------------------------------------------------------------

**Syntax**:

	dihedral_type	<type-1> <type-2> <type-3> <type-4> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines a dihedral type and its related attributes. A dihedral type is identified by
the types of the four atoms involved.

The parameters _type-1_, _type-2_, _type-3_, and _type-4_ are identifiers of previously defined atom
types. The order of such identifiers is relevant in the case of dihedral types, that is, if the same
identifiers are disposed in the inverse order, they will denote a different dihedral type. Wildcard
characters (? and *) can be used to refer to multiple atom types. If a type-related [prefix/suffix]
has been previously activated, then each actual atom type identifier will contain such prefix and/or
suffix added to _type-x_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except
comment tags (#). When a [LAMMPS data file] is generated using the command [write], then the content
of _attribute-list_ is assumed to be the parameters of the [LAMMPS dihedral style] associated with
the dihedral type in question.

It is possible to define multiple identical dihedral types. This can be done, for instance, when the
torsional potential is modeled as a sum of multiple terms, each one corresponding to a different
[LAMMPS dihedral style]. Note that two or more dihedral types can be defined identically, even with
the use of wildcard characters. However, they cannot be defined ambiguously.

Playmol detects dihedrals automatically when chemical bonds are created. If a detected dihedral fits
to a previously defined dihedral type, then Playmol will add it to a list of dihedrals to be used
later on, for instance, in the creation of a [LAMMPS data file]. If the corresponding dihedral type
does not exist at the moment a new dihedral is detected, then Playmol will ignore such dihedral and
produce a warning message. Besides the automatic detection, a dihedral can also be explicitly
created using the [extra] command.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dihedral_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the given attributes will be associated with every detected or specified
dihedral that involves four consecutive atoms whose type identifiers are _CH3_ _CH2_ _CH2_ _CH2_,
in this specific order.

**See also**:

[atom_type], [write], [prefix/suffix], [extra]

----------------------------------------------------------------------------------------------------
<a name="improper_type"></a>
improper_type
----------------------------------------------------------------------------------------------------

**Syntax**:

	improper_type	<type-1> <type-2> <type-3> <type-4> <attribute-list>

* _type-x_ = identifier of a previously defined atom type
* _attribute-list_ = list of attributes related to the bond type being defined

**Description**:

This command defines an improper type and its related attributes. An improper type is identified by
the types of the four atoms involved.

The parameters _type-1_, _type-2_, _type-3_, and _type-4_ are identifiers of previously defined atom
types. The order of such identifiers is relevant in the case of improper types, that is, if the same
identifiers are disposed in the inverse order, they will denote a different improper type. Wildcard
characters (? and *) can be used to refer to multiple atom types. If a type-related [prefix/suffix]
has been previously activated, then each actual atom type identifier will contain such prefix and/or
suffix added to _type-x_.

The parameter _attribute-list_ can be a sequence of strings containing any characters, except
comment tags (#). When a [LAMMPS data file] is generated using the command [write], then the content
of _attribute-list_ is assumed to be the parameters of the [LAMMPS improper style] associated with
the improper type in question.

It is possible to define multiple identical improper types. This can be done, for instance, when the
torsional potential is modeled as a sum of multiple terms, each one corresponding to a different
[LAMMPS improper style]. Note that two or more improper types can be defined identically, even with
the use of wildcard characters. However, they cannot be defined ambiguously.

Impropers can be manually created or detected using the command [improper], which requires a
previously defined improper type.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
improper_type	CH3 CH2 CH2 CH2 1.41095 -0.27100 3.14484 0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the given attributes will be associated with every detected or specified
improper that involves four consecutive atoms whose type identifiers are _CH3_ _CH2_ _CH2_ _CH2_,
in this specific order.

**See also**:

[atom_type], [improper], [write], [prefix/suffix]

----------------------------------------------------------------------------------------------------
<a name="atom"></a>
atom
----------------------------------------------------------------------------------------------------

**Syntax**:

	atom		<id> <type> [<charge>]

* _id_ = a name to be assigned to the atom being created
* _type_ = the name of a previously defined atom type
* _charge_ (optional) = charge of the specified atom

**Description**:

This command creates an atom of a given type.

The parameter _id_ must be a single character string with no comment tags (#) and no wildcard
characters (* or ?). If an atom-related [prefix/suffix] has been previously activated, then the
actual atom identifier will contain such prefix and/or suffix added to _id_. The same identifier
cannot be assigned to distinct atoms.

The parameter _type_ is the identifier of a previously defined atom type. A unique identifier must
be provided, with no use of wildcard characters (* or ?). If a type-related [prefix/suffix] has been
previously activated, then each actual atom type identifier will contain such prefix and/or suffix
added to _type-x_. The specified atom type must have a [mass] already define to it.

The optional parameter _charge_ is the total or partial charge of the specified atom. Note that one
can omit this parameter and assign the charge afterwards using the command [charge], which has the
advantage of admitting wildcards in order to perform multiple charge assignments simultaneously.

When an atom is created, a new monoatomic molecule is automatically detected. If there are, for
instance, _N_ already existing molecules, the new monoatomic molecule will receive an index _N+1_.
If a [bond] involving this atom is created afterwards, then such monoatomic molecule will no longer
exist.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atom	C1 C
atom	C2 C
atom    H1 H
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**See also**:

[atom_type], [mass], [bond], [charge]

----------------------------------------------------------------------------------------------------
<a name="charge"></a>
charge
----------------------------------------------------------------------------------------------------

**Syntax**:

	charge		<atom> <value>

* _atom_ = the name of an atom
* _value_ = charge of the specified atom

**Description**:

This command defines the electric charge (total or partial) of an atom.

The parameter _atom_ is an [atom] identifier. If an atom-related [prefix/suffix] has been previously
activated, then the actual atom identifier will contain such prefix and/or suffix added to _atom_.
Wildcard characters (* or ?) can be used to simultaneously assign the same charge to multiple atoms.
__Important__: the presence of wildcards makes a charge value applicable to atoms defined either
beforehand or afterwards.

If an [atom] has no explicitly defined charge, then it will be considered as neutral.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
charge		O* -1.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, a change value of _-1.0_ is assigned to all atoms whose identifiers start with
the letter _O_.

**See also**:

[atom]

----------------------------------------------------------------------------------------------------
<a name="bond"></a>
bond
----------------------------------------------------------------------------------------------------

**Syntax**:

	bond		<atom-1> <atom-2> [<atom-3> <atom-4> ...]

* _atom-x_ = the name of a previously defined atom

**Description**:

This command creates chemical bonds between atoms.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be
provided, with no use of wildcard characters (* or ?).

Parameters _atom-1_ and _atom-2_ are mandatory. Additional parameters _atom-x_ (for _x_ > 2) are
optional. In any case, _atom-1_ is the central atom to which all other atoms will be connected. If
an atom-related [prefix/suffix] has been previously activated, then the actual atom identifier will
contain such prefix and/or suffix added to _atom-x_.

If a bond is created between atoms that belong to the same molecule, then a cyclic structure is
formed and the number of molecules does not change. However, if the two atoms belong to distinct
molecules, then these molecules are fused together and the number of molecules is decreased.
Supposing that the fused molecules have indices _i_ and _j_, the index of the new molecule will be
_min(i,j)_. In order to keep a consistent sequence, all existing molecules whose indexes are greater
than _max(i,j)_ will have their indexes decreased in one unit.

New angles and new dihedrals are detected automatically when a bond is created. If a detected
structure fits to a previously defined type, then Playmol will add it to a list to be used later
on, for instance, in the creation of a [LAMMPS data file]. If the corresponding type does not exist
at the moment a new structure is detected, then Playmol will ignore such structure and produce a
warning message.

After chemical bonds have been created, one can visualize the current list of molecules by invoking
the command [write] with the _summary_ format.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bond	O H1 H2
write	summary
quit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, two chemical bonds _O_-_H1_ and _O_-_H2_ are created. Then, a command [write]
is issued and Playmol execution is interrupted, so that the user can visualize the updated list of
molecules.

**See also**:

[bond_type], [atom_type], [write], [quit]

----------------------------------------------------------------------------------------------------
<a name="improper"></a>
improper
----------------------------------------------------------------------------------------------------

**Syntax**:

	improper	<atom-1> <atom-2> <atom-3> <atom-4>

or

	improper	search

* _atom-x_ = name of a previously defined atom

**Description**:

This command creates an improper involving the specified atoms or search for impropers.

__Important__: Unlike angles and proper dihedrals, impropers are not detected automatically when the
command [bond] is executed. This is so because improper dihedrals can be defined in several distinct
ways.

Each parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be
provided, with no use of wildcard characters (* or ?). If an atom-related [prefix/suffix] has been
previously activated, then the actual atom identifier will contain such prefix and/or suffix added
to _atom-x_.

When using explicit improper definition (first syntax above), a previously defined [improper_type]
matching the four atoms in the specified order is required.

Using the keyword _search_ (second syntax above), Playmol will search for impropers composed of any
four atoms I, J, K, and L in which atom K is simultaneously bonded to atoms I, J, and L. A detected
improper will only be effectively created if the corresponding [improper_type] has been previously
defined.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
improper_type	HC HC C HC	1.41095 -0.27100 3.14484 0
improper     	search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, an improper type is defined for all-atom methyl groups and then Playmol
searches for impropers of this type.

**See also**:

[improper_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="extra"></a>
extra
----------------------------------------------------------------------------------------------------

**Syntax**:

	extra	bond		<atom-1> <atom-2>

or

	extra	angle		<atom-1> <atom-2> <atom-3>

or

	extra	dihedral	<atom-1> <atom-2> <atom-3> <atom-4>

* _atom-x_ = name of a previously defined atom

**Description**:

This command creates an extra structure (bond, angle, dihedral) involving the specified atoms, which
must all belong to the same molecule.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be
provided, with no use of wildcard characters (* or ?). A corresponding structure type (that is,
[bond_type], [angle_type], or [dihedral_type]) must have been previously defined. If an atom-related
[prefix/suffix] has been previously activated, then the actual atom identifier will contain such
prefix and/or suffix added to _atom-x_.

__Important__: Playmol does not automatically search for angles and dihedrals when a chemical bond
is created with the command [extra].

NOTE: extra structures are only required by some molecular models that define unconventional bonds,
angles, or dihedrals, such as those that involve atoms which are not actually bonded in a molecule,
for example. Most models, however, will involve only automatically detectable angles and dihedrals.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extra bond	C1 C2
extra angle	C1 C2 C3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**See also**:

[bond_type], [angle_type], [dihedral_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="link"></a>
link
----------------------------------------------------------------------------------------------------

**Syntax**:

	link	<atom-1> <atom-2>

* _atom-x_ = name of a previously defined atom

**Description**:

This command creates a virtual link between the specified atoms. If these atoms belong to distinct
molecules, then such molecules will be fused together, but without actually creating a chemical bond
between the linked atoms.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be
provided, with no use of wildcard characters (* or ?). If an atom-related [prefix/suffix] has been
previously activated, then the actual atom identifier will contain such prefix and/or suffix added
to _atom-x_.

NOTE: sometimes, it is useful to consider two molecules as if they were a single one. For instance,
when ions are present, one might want to consider an anion/cation pair as a single rigid structure,
so that they can be packed together using the [packmol] command. After that, one can remove the
link using the [unlink] command before writing down the produced configuration as, for instance, a
[LAMMPS data file].

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
link	Na+ Cl-
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**See also**:

[atom], [packmol], [unlink]

----------------------------------------------------------------------------------------------------
<a name="unlink"></a>
unlink
----------------------------------------------------------------------------------------------------

**Syntax**:

	unlink	<atom-1> <atom-2>

* _atom-x_ = name of a previously defined atom

**Description**:

This command removes a previously defined link between the specified atoms. If such link is the only
connection between two otherwise independent parts of a molecule, then such molecule will be
disassembled accordingly. If necessary, the list of atomic coordinates of instantiated molecules
will be rearranged in order to keep contiguousness for the atoms of each molecule.

The parameter _atom-x_ is the identifier of a previously created atom. A unique identifier must be
provided, with no use of wildcard characters (* or ?). If an atom-related [prefix/suffix] has been
previously activated, then the actual atom identifier will contain such prefix and/or suffix added
to _atom-x_. The two identifiers must correspond to a previously created [link].

NOTE: sometimes, it is useful to consider two molecules as if they were a single one. For instance,
when ions are present, one might want to consider an anion/cation pair as a single rigid structure,
so that they can be packed together using the [packmol] command. After that, one can remove the
link using the [unlink] command before writing down the produced configuration as, for instance, a
[LAMMPS data file].

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unlink	Na+ Cl-
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**See also**:

[atom], [packmol], [link]

----------------------------------------------------------------------------------------------------
<a name="build"></a>
build
----------------------------------------------------------------------------------------------------

**Syntax**:

	build     [<file>]

* _file_ (optional) = name of a file containing geometric information

**Description**:

This command reads geometric information and uses them to determine atomic coordinates for one or
more molecules.

The optional parameter _file_ is the name of a text file containing geometric information, whose
required format is described below. If it is omitted, then the geometric information must be
provided in the subsequent lines. In this sense, issuing the command `build <file>` is equivalent to
issuing `build`, followed by `include <file>`.

The geometric information that Playmol expects to receive has the following format:

* The first non-empty line must contain the number _N_ of geometric assignments that will be made.
* The subsequent non-empty lines must contain the _N_ geometric assignments, which can possibly be
intertwined with [define] commands, [include] commands, [for/next] constructs, and [if/then/else]
constructs. No other commands are permitted until all _N_ assignments have been made. Each geometric
assignments consists of an [atom] identifier followed by some fields separated by spaces and/or
tabs. Only the following formats are valid:

Coordinates _x_, _y_, and _z_ of _atom-I_:

	<atom-I>  <x>  <y>  <z>

The _length_ of a bond formed by atoms _I_ and _J_:

	<atom-I>  <atom-J>  <length>

The _length_ of a bond formed by atoms _I_ and _J_ and the _angle_ (in degrees) between the bonds
formed by atoms _I_, _J_, and _K_:

	<atom-I>  <atom-J>  <length>  <atom-K>  <angle>

The _length_ of a bond formed by atoms _I_ and _J_, the _angle_ (in degrees) between the bonds
formed by atoms _I_, _J_, and _K_, and the _torsion_ angle (in degrees) of the proper dihedral
formed by atoms _I_, _J_, _K_, and _L_:

	<atom-I>  <atom-J>  <length>  <atom-K>  <angle>  <atom-L>  <torsion>

An active atom-related [prefix/suffix] will apply to every parameter _atom-x_ present in the
assignments described above. In all cases, atom _I_ is the one whose coordinates are being assigned,
while atoms _J_, _K_, and _L_ must be atoms whose coordinates have already been assigned (please see
the examples below).

Every molecule must be instantiated at once, meaning that the geometric information regarding all
atoms of a given molecule must be provided contiguously. However, no particular order is required.
Multiple molecules of one or more compounds can be instantiated either using separate _build_
constructs or a unique _build_ constructs with all atoms together.

Note that the ubiquitous [xyz file format] can be used to provide geometric information for the
_build_ constraint, given that:

* The first line contains the number of atomic coordinates (_N_), as usual.
* The second line is empty or starts with a comment mark (#).
* The _N_ subsequent lines contain the following four fields separated by spaces and/or tabs:

	<atom-id>  <x>  <y>  <z>

Therefore, one can use an external software to generate the atomic coordinates (e.g., [Avogadro]),
save them in the [xyz file format], and replace the usual element symbols by atom identifiers.

The list of molecular structures created via _build_ can be employed later on as prototypes for
replication and packing with the [packmol] command. They can also be directly employed to produce a
simulation box when the command [write] is invoked. In this case, the origin of the Cartesian space
will be located at the center of the simulation box.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build	H2O.info
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the geometric information of one or more copies of a water molecule is read
from the file _H2O.info_.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build
3
O   0.00000  0.00000  0.00000
H1  O  0.9572
H2  O  H1  0.9572  104.52
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, the atomic coordinates of a water molecule are guessed by positioning the
oxygen atom at the origin, defining the lengths of the oxygen-hydrogen bonds as 0.9572 Å, and then
defining the hydrogen-oxygen-hydrogen angle as 104.52 degrees.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build
5
C1              0.0     0.0    0.0
C2 C1           1.54
C3 C2 C1        1.33    114
C4 C3 C2 C1     1.54    114      0
C5 C4 C3 C2     1.54    114    180

define	file as cis_2_pentene.lammpstrj
box   	density 1.0
write 	lammpstrj ${file}
shell 	gnome-terminal -e "vmd ${file}"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example above, a united-atom cis-2-pentene molecule (C1-C2=C3-C4-C5) is defined. The length
of the double bond is 1.33 Å and that of the single bonds is 1.54 Å. All bond-bending angles are 114
degrees. The torsional angle of the dihedral C1-C2=C3-C4 is zero (cis) and that of the dihedral
C2=C3-C4-C5 is 180 degrees (trans). In the end, a configuration file is written down and [vmd] is
involved in a new terminal window for visualization.

**See also**:

[atom], [bond], [box], [write], [include]

----------------------------------------------------------------------------------------------------
<a name="prefix_suffix"></a>
prefix/suffix
----------------------------------------------------------------------------------------------------

**Syntax**:

	prefix		<target> <string>

or

	suffix		<target> <string>

* _target_ = _types_ or _atoms_
* _string_ = _none_ or a character string to be used as prefix or suffix for atom types or atoms

**Description**:

This commands activates/deactivates a prefix or a suffix to be added to every [atom_type] identifier
or to every [atom] identifier in subsequent commands.

The parameter _target_ indicates which type of identifier will be modified by the activated prefix
or suffix. If _target_ = _types_, then the prefix/suffix will be applied to every [atom_type]
identifier in subsequent commands. Otherwise, if _target_ = _atoms_, then the prefix/suffix will be
applied to every [atom] identifier in subsequent commands.

The parameter _string_ must be a single character string with no comment tags (#) and no wildcard
characters (* or ?). Such string will be used in subsequent commands as a prefix or suffix for
[atom_type] identifiers or for [atom] identifiers, depending on the specified parameter _target_.
If _string_ = _none_, then no prefix/suffix will be added for the specified _target_ in subsequent
commands. Therefore, using _string_ = _none_ deactivates a currently active prefix or suffix.

Prefixes and suffixes can be simultaneously active.

The following commands are affected by the `prefix types` and `suffix types` commands:

* [atom_type], [mass], [diameter], [bond_type], [angle_type], [dihedral_type], [improper_type], and
[atom].

The following commands are affected by the `prefix atoms` and `suffix atoms` commands:

* [atom], [bond], [charge], [improper], [extra], [link], and [build].

__Note__: an active, atom-related prefix/suffix will NOT be automatically applied to arguments of
the `mol()` function. Therefore, molecule specifications within the commands [align] and [packmol]
require manual inclusion of prefixes and/or suffixes.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prefix		atoms	H2O-
suffix		types	A
atom		H1 H
atom		H2 H
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above creates two atoms whose identifiers are _H2O-H1_ and _H2O-H2_. Both atoms are of
the same type _HA_.

**See also**:

[atom_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="box"></a>
box
----------------------------------------------------------------------------------------------------

**Syntax**:

	box		lengths <lx> <ly> <lz>

or

	box		angles <alpha> <beta> <gamma>

or

	box		volume <value> [aspect <ax> <ay> <az>]

or

	box		density <value> [aspect <ax> <ay> <az>]

* _lx_, _ly_, _lz_ = lengths of the box edges in directions x, y, and z
* _alpha_, _beta_, _gamma_ = angles between the box edges
* _value_ = density of the simulation box
* _ax_, _ay_, _az_ (optional) = Scaling factors in directions x, y, and z

**Description**:

This command specifies the side lengths, the angles, the volume, or the density of a simulation box.

The parameters _lx_, _ly_, and _lz_ are the lengths of the simulation box edges in the three spatial
directions.

The parameters _alpha_, _beta_, and _gamma_ are the angles, in degrees, between edge vectors of a
non-orthogonal box (_alpha_ between _ly_ and _lz_, _beta_ between _lx_ and _lz_, and _gamma_ between
_lx_ and _ly_). The command `box angles` requires that a command `box lengths` has been previously
issued.

The parameter _value_ is either the volume or the density of the simulation box, depending on the
kind of specification.

The optional keyword _aspect_ with parameters _ax_, _ay_, and _az_ can be used to determine the
relative box sizes in the three spatial directions (the absolute values of these three parameters
are meaningless). If the keyword _aspect_ is omitted, Playmol will consider a cubic box.

__Important__: when the box density is specified, Playmol will automatically calculate the box side
lengths when necessary. For this, it will sweep the whole list of molecular structures defined via
[build] and [packmol] commands in order to determine the total mass of the system. Playmol will not
perform any unit conversion. Consequently, the specified density value must be consistent with the
mass and length units used in other commands. For instance, if [LAMMPS real units] are considered
(mass in g/mol and lengths in Å), then the density must be provided in Da/Å³ (Daltons per cubic
Angstrom). In this case, values originally in g/cm³ must be multiplied by 0.60221413.

When building a simulation box, Playmol considers that its geometric center is in the origin of the
Cartesian space. Playmol does not check whether the previously specified atomic coordinates fit
inside this box.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box	lengths 30 30 30
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above specifies a cubic simulation box.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box	volume 27000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above is completely equivalent to the preceding one, since the default aspect is cubic.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box	density 0.602214 aspect 1 1 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above specifies a simulation box with density equal to 0.602214 in [LAMMPS real units]
(equivalent to 1.0 g/cm³). The box is twice as elongated in the z direction as in the other two
directions.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box	lengths 28.224 28.224 34.443
box	angles	90 90 120
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above specifies a non-orthogonal simulation box that can be used to simulate a hexagonal
lattice system (_lx_ = _ly_ ≠ _lz_ and _alpha_ = _beta_ ≠ _gamma_).

**See also**:

[build], [packmol], [write]

----------------------------------------------------------------------------------------------------
<a name="align"></a>
align
----------------------------------------------------------------------------------------------------

**Syntax**:

	align		 <molecule> <axis-1> <axis-2>

* _molecule_ = the index of an existing molecule with predefined coordinates
* _axis-1_ = _x_, _y_, or _z_
* _axis-2_ = _x_, _y_, or _z_ (different from _axis-1_)

**Description**:

This command aligns the principal axes of a given molecule with the system's Cartesian axes. The
principal axes are defined as those for which the molecule's inertia tensor is diagonal.

The parameter _molecule_ specifies the molecule to be aligned with the Cartesian axes. This
specification can be done in either of the following ways:

1. By directly employing the numerical index of the compound. The command [write], with its option
_summary_, can be helpful for checking the indexes of the existing ones.

2. By using the function `mol(atom)`, where _atom_ is the identifier of an existing [atom]. Such
identifier must be tightly placed inside the parentheses (that is, without any spaces and/or tabs).
Note that an active, atom-related [prefix/suffix] will NOT be applied automatically to _atom_.

Except in a few special cases, the specified molecule must have already been instantiated (see the
[Playmol Basics] section). In this case, the first set of atomic coordinates previously provided
for such molecule will be employed. The special cases, which might not require previous
instantiation, are: (1) monoatomic molecules, (2) diatomic molecules whose involved [bond_type]
describes a harmonic potential, and (3) triatomic molecules whose involved [bond_type] and
[angle_type] both describe harmonic potentials. In these cases, Playmol will consider the second
[bond_type]'s attribute as an equilibrium distance and the second [angle_type]'s attribute as an
equilibrium angle (in degrees).

The parameter _axis-1_ must be equal to _x_, _y_, or _z_, and specifies the Cartesian axis with
which the molecule will be aligned considering its most elongated principal axis.

The parameter _axis-2_ must also be equal to _x_, _y_, or _z_, but different from _axis-1_. It
specifies the Cartesian axis with which the molecule will be aligned considering its second most
elongated principal axis. Of course, the shortest principal axis will be aligned to the Cartesian
axis not explicitly specified in the command _align_.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
align		1 x z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above aligns the most elongated principal axis of molecule 1 to the Cartesian axis _x_,
and its second most elongated principal axis to axis _z_.

**See also**:

[write], [packmol]

----------------------------------------------------------------------------------------------------
<a name="packmol"></a>
packmol
----------------------------------------------------------------------------------------------------

**Syntax**:

	packmol		<keyword> <arguments> [<keyword> <arguments> ...]

* _keyword_ =  _seed_ or _diameter_ or _nloops_ or _retry_ or _fix_ or _copy_ or _pack_ or _action_

**Keywords and arguments**:

* seed <iseed>
    * _iseed_ = an integer seed for Packmol's random number generator (_default_ = 1234)
* diameter <D>
    * _D_ = diameter assigned to all atoms expect those with explicit definition (_default_ = 2.5)
* nloops <N>
    * _N_ = the maximum number of iterations of the Packmol algorithm (_default_ = 50)
* retry <factor>
    * _factor_ = a scaling factor for the tolerance in successive packing attempts (_default_ = 1.0)
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

This command invokes the [Packmol package] to:

1. Build an overlap-free molecular packing by considering atoms as hard spheres. In this case, an
overlap is considered to occur if the distance between any two atoms of distinct molecules is
smaller than the sum of their radii. After a packing is successfully built, all atomic coordinates
defined beforehand are deleted (as if a [reset] command was issued with option _xyz_) and fully
replaced with the new Packmol-generated coordinates.

2. Create input files so that the molecular packing can be produced by executing Packmol externally.

The keywords _seed_, _diameter_, and _nloops_ change the values of some parameters that affect the
behavior of Packmol.

The parameter _iseed_ is a seed for Packmol's pseudo-random number generator. Even with all other
parameters unchanged, different values of _iseed_ will produce distinct random packings. It must be
an integer number and its default value is _1234_.

The parameter _D_ is the default atom diameter, that is, the diameter to be assigned to all atom
types except those whose diameters have been explicitly defined using the [diameter] command. It
must be a nonnegative real number and its default value is _2.5_. This might be a suitable value for
all-atom models when distances are expressed in Angstroms.

The parameter _N_ is the maximum number of iterations allowed for the Packmol algorithm. If this
number is exceeded, then a packing attempt is considered as unsuccessful. It must be a positive
integer number and its default value is _50_.

The parameter _factor_ is a diameter reduction factor. It must satisfy 0.0 &lt; _factor_ &le; 1.0
and its default value is 1.0. If _factor_ = 1.0, then Packmol will be invoked only once with the
specified atom diameters and will produce an error message if the packing attempt fails. Otherwise,
Packmol will be invoked as many times as necessary to achieve a successful packing, with all atom
diameters being iteratively multiplied by the value of _factor_ at each new attempt.

The keywords _fix_, _copy_, and _pack_ require specifying an existing molecule. Except in a few
special cases, such molecule must have already been instantiated (see the [Playmol Basics] section).
In this case, the first sets of atomic coordinates previously provided for such molecule will be
used to define a rigid-body that can undergo replications, translations, and rotations. The special
cases, which might not require previous instantiation, are: (1) monoatomic molecules, (2) diatomic
molecules whose involved [bond_type] describes a harmonic potential, and (3) triatomic molecules
whose involved [bond_type] and [angle_type] both describe harmonic potentials. In these cases,
Playmol will consider the second [bond_type]'s attribute as an equilibrium distance and the second
[angle_type]'s attribute as an equilibrium angle (in degrees).

A molecule specification for the _fix/copy/pack_ keywords can be done in either of the following
ways:

1. By directly employing the numerical index of the compound. The command [write], with its option
_summary_, can be helpful for checking the indexes of the existing ones. One may note, however,
that these indexes can possibly change during the execution of a script (please see the [bond]
command description). Thus, the specification must refer to the index which the desired compound
will have at the moment of a `packmol action` command (see below).

2. By using the function `mol(atom)`, where _atom_ is the identifier of an existing [atom]. Such
identifier must be tightly placed inside the parentheses (that is, without any spaces and/or tabs).
At the moment of a `packmol action` command, the function will return the index of the compound that
contains the specified atom. Note that an active, atom-related [prefix/suffix] will NOT be applied
automatically to _atom_.

The usage of each keyword _fix_, _copy_, or _pack_ is:

	fix <molecule> <x> <y> <z>

This makes a copy of a molecule and, while keeping the original orientation, places its geometric
center at the provided coordinate (_x_, _y_, _z_).

	copy <molecule> <N>

This makes _N_ copies of a molecule in random positions, but keeping the original orientation
(i.e. all _N_ copies will be aligned to each other in the final packing). This is useful for
facilitating the packing of elongated molecules.

	pack <molecule> <N>

This makes _N_ copies of a molecule in random positions and with random orientations.

The keyword _action_ is used to create Packmol input files or to execute Packmol. The following
options are available:

	action execute

This option executes Packmol in order to build the desired molecular packing. It requires the
previous definition of a simulation [box]. Moreover, it requires that at least one keyword _fix_
_copy_, or _pack_ has appeared in a previous packmol command or appears in the same packmol command,
either before or after the keyword _action_. If the parameter _retry_ is currently equal to 1.0 (its
default value), then Packmol will do only one packing attempt with the specified _tolerance_ and
produce a warning message in case such attempt fails. If _retry_ is smaller than 1.0, then Packmol
will keep trying until a successful attempt is achieved, with _tolerance_ being iteratively
multiplied by the _retry_ value at each new attempt. IMPORTANT: if the packing succeeds, then the
current list of atomic coordinates is replaced by the new coordinates generated by Packmol. After
that, the command [write] can be used to create a LAMMPS configuration file with the attained
packing.

	action setup

This option generates an input file named _packmol.inp_ and coordinate files _molecule-x.inp_, where
_x_ is an index for each involved molecule. These files are prepared for running Packmol externally
in order to generate a file _packmol-output.xyz_ containing the final packing if the algorithm
succeeds with the given tolerance. For illustration, one may notice that a successful use of the
packmol command with options _retry 1.0 action execute_ would have exactly the same result as the
following sequence of commands:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
packmol 	action setup
shell   	packmol < packmol.inp
reset		xyz
build     	packmol_output.xyz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nevertheless, the real usefulness of the option _setup_ is to allow the file _packmol.inp_ to be
manually edited so that some additional constraints can be imposed. The user is referred to
[Packmol User's Guide] for additional information.

The Packmol algorithm is described in the following paper:

    @article{packmol_paper,
        author = {Martinez, L. and Andrade, R. and Birgin, E. G. and Martinez, J. M.},
        title = {PACKMOL: A package for building initial configurations for molecular dynamics
                 simulations},
        journal = {J. Comput. Chem.},
        volume = {30},
        number = {13},
        pages = {2157--2164},
        year = {2009}
        doi = {10.1002/jcc.21224},
    }

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box    		density 0.602214
packmol		tolerance 3.0 retry 0.9 fix 1 0.0 0.0 0.0 pack 2 1000
packmol		action execute
write		lammps system.data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above uses Packmol to create a random packing of molecules with density equal to
0.602214 Da/Å³ (1.0 g/cm³) in which one copy of molecule 1 is centered at the origin and 1000 copies
of molecule 2 are packed with random positions and random orientations. The desired mininum
intermolecular atomic distance is initially set to 3.0 Å, but Packmol will keep reducing the
tolerance in 90% until the packing is successful. Finally, a [LAMMPS data file] is generated.

**See also**:

[reset], [diameter], [box], [build], [write], [bond_type], [angle_type]

----------------------------------------------------------------------------------------------------
<a name="write"></a>
write
----------------------------------------------------------------------------------------------------

**Syntax**:

	write		 <format> [<file>]

* _format_ = _playmol_ or _lammps_ or _summary_ or _xyz_ or _lammpstrj_
* _file_ (optional) = name of a file to be created

**Description**:

This command writes down the molecular system in a file or in the standard output device.

The parameter _format_ must be one of the following options:

* __playmol__: the output will contain Playmol commands that could be used in another script to
build the same system. For illustration, detected angles and dihedrals appear as commented lines.
Type and atom prefixes are explicitly added to the corresponding identifiers.

* __lammps__: the command will produce information in the LAMMPS configuration file format, which
can be used as an initial configuration for a Molecular Dynamics simulation using LAMMPS through its
command [read_data].

* __summary__: this option will print a summary of the system characteristics, including the amount
of every defined and detected structure such as angles, dihedrals, and molecules. Properties of each
molecular species will also be printed, such as its mass, its atoms, and the number of defined sets
of coordinates. This is useful for debugging purposes.

* __xyz__: writes down the list of atomic coordinates using the [xyz file format], but with element
symbols replaced by atom identifiers. This is useful for using with another Playmol script or for
visualization purposes.

* __lammpstrj__: writes down the list of atomic coordinates using the LAMMPS trajectory format. This
is useful for visualization with [VMD].

The optional parameter _file_ is the name of the file which will contain the system info. If it is
omitted, the info will be written in the standard output unit (the computer screen, in most cases).

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write	playmol water.mol
write	lammps water.data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

----------------------------------------------------------------------------------------------------
<a name="include"></a>
include
----------------------------------------------------------------------------------------------------

**Syntax**:

	include		<file>

* _file_ = name of a file containing additional playmol commands

**Description**:

This command redirects Playmol to read and execute commands from another script.

**Examples**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include		H2O.mol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**See also**:

[atom_type], [atom]

----------------------------------------------------------------------------------------------------
<a name="reset"></a>
reset
----------------------------------------------------------------------------------------------------

**Syntax**:

	reset		<list>

* _list_ = *all* or *atom* or *charge* or *bond* or *link* or *angle* or *dihedral* or *improper* or
*extra_bond* or *extra_angle* or *extra_dihedral* or *xyz* or *packmol*

**Description**:

This command resets one or more lists of entities previously created. The resettable lists are those
in the table below.

| Index | List                                 |
|------:|:-------------------------------------|
|     1 | atom types                           |
|     2 | masses                               |
|     3 | diameters                            |
|     4 | bond types                           |
|     5 | angle types                          |
|     6 | dihedral types                       |
|     7 | improper types                       |
|     8 | atoms                                |
|     9 | charges                              |
|    10 | bonds                                |
|    11 | angles                               |
|    12 | dihedrals                            |
|    13 | extra bonds                          |
|    14 | extra angles                         |
|    15 | extra dihedrals                      |
|    16 | impropers                            |
|    17 | links                                |
|    18 | atomic coordinates                   |
|    19 | packmol fix/copy/pack specifications |

The parameter _list_ must be one of the following options:

* __all__: deletes all lists previously created (_1_ to _19_).
* __atom__: deletes the list of atoms and all its dependent lists (_8_ to _19_).
* __charge__: deletes the list of atomic charges (_9_).
* __bond__: deletes the lists from _10_ to _17_ in the table above and resets all atoms as
monoatomic molecules. The lists of coordinates and packmol specifications remain unchanged.
* __angle__: deletes the lists of angles and extra angles (_11_ and _14_).
* __dihedral__: deletes the lists of dihedrals and extra dihedrals (_12_ and _15_).
* __extra_bond__: deletes the list of extra bonds (_13_).
* __extra_angle__: deletes the list of extra angles (_14_).
* __extra_dihedral__: deletes the list of extra dihedrals (_15_).
* __improper__: deletes the list of impropers (_16_).
* __link__: deletes the list of virtual links (_17_) and splits all molecules accordingly, as if
the [unlink] command was issued individually for every link present in the list.
* __xyz__: deletes the list of atomic coordinates (_18_).
* __packmol__: deletes the list of [packmol fix/copy/pack specifications (_19_).

**See also**:

[packmol], [unlink]

----------------------------------------------------------------------------------------------------
<a name="shell"></a>
shell
----------------------------------------------------------------------------------------------------

**Syntax**:

	shell		<command>

* _command_ = an external shell command

**Description**:

This command executes an external shell command.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shell	mv packmol.inp system.inp
shell	packmol < system.inp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above uses Linux command _mv_ to rename a file from _packmol.inp_ to _system.inp_ and
then executes Packmol using it as input.

**See also**:

[packmol]

----------------------------------------------------------------------------------------------------
<a name="quit"></a>
quit
----------------------------------------------------------------------------------------------------

**Syntax**:

	quit [all]

**Description**:

This command interrupts the execution of a Playmol script, thus being useful for debugging purposes.

The optional keyword _all_ forces Playmol to stop completely. In the absence of the keyword, the
command _quit_ has the same behavior of a normal end of file. For instance, if the script has been
invoked by another script via the [include] command, execution of the invoking script continues.

**Example**:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write summary
quit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example above writes a summary of the current molecular system and then quits Playmol.

**See also**:

[include]

----------------------------------------------------------------------------------------------------

[define]:		#define
[for/next]:		#for_next
[if/then/else]:		#if_then_else
[atom_type]:		#atom_type
[mass]:			#mass
[diameter]:		#diameter
[bond_type]:		#bond_type
[angle_type]:		#angle_type
[dihedral_type]:	#dihedral_type
[improper_type]:	#improper_type
[atom]:			#atom
[charge]:		#charge
[bond]:			#bond
[extra]:		#extra
[link]:			#link
[unlink]:		#unlink
[improper]:		#improper
[build]:		#build
[box]:			#box
[packmol]:		#packmol
[align]:		#align
[write]:		#write
[prefix/suffix]:	#prefix_suffix
[include]:		#include
[reset]:		#reset
[shell]:		#shell
[quit]:			#quit

[Playmol Basics]:		basics.html
[LAMMPS real units]:		http://lammps.sandia.gov/doc/units.html
[LAMMPS data file]:		http://lammps.sandia.gov/doc/read_data.html
[LAMMPS pair style]:		http://lammps.sandia.gov/doc/pair_style.html
[LAMMPS bond style]:		http://lammps.sandia.gov/doc/bond_style.html
[LAMMPS angle style]:		http://lammps.sandia.gov/doc/angle_style.html
[LAMMPS dihedral style]:	http://lammps.sandia.gov/doc/dihedral_style.html
[LAMMPS improper style]:	http://lammps.sandia.gov/doc/improper_style.html
[read_data]:			http://lammps.sandia.gov/doc/read_data.html
[xyz file format]:		http://openbabel.org/wiki/XYZ_(format)
[Packmol package]:		http://www.ime.unicamp.br/~martinez/packmol
[Packmol User's Guide]:		http://www.ime.unicamp.br/~martinez/packmol/quickguide
[VMD]:				http://www.ks.uiuc.edu/Research/vmd
[Avogadro]:			http://avogadro.cc/wiki/Main_Page
[ASCII]:			https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart
