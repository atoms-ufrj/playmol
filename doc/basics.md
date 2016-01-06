Playmol Basics      {#basics}
==============

Introduction
----------------------------------------------------------------------------------------------------

Playmol is a molecular model builder which intends to facilitate the construction of initial setups
for Molecular Dynamics simulations of dense systems. It was initially designed as a bridge between
the molecular packing generator [Packmol] and the popular Molecular Dynamics simulator [LAMMPS].

The principle of Playmol consists in the bottom-up assembly of molecules and molecular systems by a
series of steps like:

1. __Classification__: atom types are defined and combinations of these types are registered as
templates for the identification of chemical bonds, angles, dihedrals, and impropers. These are the
topological features most usually considered in classical force fields, whose parameters can also be
fed to Playmol in this first step.

2. __Chemical topology definition__: atoms are created and explicitly connected via chemical bonds.
This is all Playmol needs to identify chemical compounds and automatically determine the angles and
proper dihedrals present in their molecules. Playmol does not automatically detect impropers, but
one can point them out manually or ask Playmol to search for impropers of some specific kind.

3. __Instantiation__: chemical compounds identified in the previous step remain as abstract entities
until actual molecules are positioned in space. The [build] command makes it possible to instantiate
multiple molecules of each compound in order to create a molecular system.

4. __Replication and packing__: after instantiating at least one molecule of each compound, Playmol
can invoke [Packmol] to instantiate replicas in random positions, with possibly random orientations,
while trying to avoid undesirable overlaps. This allows a dense packing of molecules to be a viable
initial configuration for molecular dynamics.

5. __Formatting and storage__: finally, the system configuration and the force field parameters are
stored in a file whose format can be understood by a simulation software such as [LAMMPS].

Practical examples of the use of Playmol can be found in the section that describes its [commands]
and in a special section with sample [scripts].

----------------------------------------------------------------------------------------------------
Playmol Scripts
----------------------------------------------------------------------------------------------------

Playmol interprets scripts composed of specific types of [commands]. Each command is a sequence of
keywords and parameter values separated by spaces and/or tabs. A script can include comments,
identified by the comment mark "#", which instructs Playmol to ignore all trailing characters in the
same line. A single command can span several lines by means of continuation marks "&", as in the
example below for the [bond] command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bond	C2       & # This is the central carbon atom,
    	C1  C3   & # which makes bonds with two other carbons
    	H21 H22    # and with two hydrogens.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The execution flow of a Playmol script is sequential, except when a [for/next] construct or an
[if/then/else] construct is found. The execution of every command will always be done in three
steps:

1. __Substitution of variables__: 

2. __Evaluation of mathematical expressions__: 

3. __Parsing and execution__: 

Note that Playmol does not perform any previous syntax checks, which means that syntax errors will
be handled at run-time.

----------------------------------------------------------------------------------------------------
<a name="variables"/>
Variables
----------------------------------------------------------------------------------------------------

Variable assignment is done via the [define] command:


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define compound as H2O
write  xyz ${compound}.xyz
diameter a b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

----------------------------------------------------------------------------------------------------
<a name="math"/>
Mathematical Expressions
----------------------------------------------------------------------------------------------------

In a Playmol script, mathematical expressions are placed between curly brackets "{" and "}". Before executing
a command line, Playmol first substitutes the value of every variable in it (please see the section
[Variables] for details), then it searches for mathematical expressions to evaluate.

<a name="Table_1"/> **Table 1**: Precedence order of mathematical operators

| Operator                                 | Description                             | Prececence  |
|:----------------------------------------:|:----------------------------------------|:-----------:|
| \|                                       | Logical OR                              | 1 (lowest)  |
| \&                                       | Logical AND                             | 2           |
| < <br> <= <br> > <br> >= <br> == <br> != | Relational operators                    | 3           |
| + <br> -                                 | Addition<br>Subtraction                 | 4           |
| * <br> / <br> %                          | Multiplication<br>Division<br>Remainder | 5           |
| ^                                        | Exponentiation                          | 6           |
| `function( )`                            | Mathematical functions ([Table 2])      | 7           |
| ( )                                      | Parentheses                             | 8 (highest) |

__Note__: There is not a specific operator for logical NOT. Its effect is obtained by using the
function `not()`.

<a name="Table_2"/> **Table 2**: Available mathematical functions

| Function | Result type      | Description                |
|:--------:|:----------------:|:---------------------------|
| abs      | Same as argument | Absolute value             |
| acos     | Real             | Inverse cosine             |
| asin     | Real             | Inverse sine               |
| atan     | Real             | Inverse tangent            |
| ceil     | Integer          | Smallest following integer |
| cos      | Real             | Cosine                     |
| cosh     | Real             | Hyperbolic cosine          |
| exp      | Real             | Exponential                |
| floor    | Integer          | Largest preceding integer  |
| int      | Integer          | Integer part               |
| ln       | Real             | Natural logarithm          |
| log      | Real             | Decimal logarithm          |
| nint     | Integer          | Nearest integer            |
| not      | 0 or 1           | Logical NOT                |
| sin      | Real             | Sine                       |
| sinh     | Real             | Hyperbolic sine            |
| sqrt     | Real             | Square root                |
| tan      | Real             | Tangent                    |
| tanh     | Real             | Hyperbolic tangent         |

----------------------------------------------------------------------------------------------------

<!-- Internal links -->
[Table 1]:		#Table_1
[Table 2]:		#Table_2
[Variables]:		#variables

<!-- External links -->
[scripts]:		scripts.html
[commands]:             commands.html
[define]:		commands.html#define
[for/next]:		commands.html#for_next
[if/then/else]:		commands.html#if_then_else
[bond]:			commands.html#bond
[build]:		commands.html#build

[LAMMPS]:		http://lammps.sandia.gov
[Packmol]:		http://www.ime.unicamp.br/~martinez/packmol
[read_data]:		http://lammps.sandia.gov/doc/read_data.html
[xyz file format]:	http://openbabel.org/wiki/XYZ_(format)
[Packmol User's Guide]:	http://www.ime.unicamp.br/~martinez/packmol/quickguide/
[VMD]:			http://www.ks.uiuc.edu/Research/vmd/

