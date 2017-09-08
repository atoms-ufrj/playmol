Playmol Basics      {#basics}
==============

Introduction
----------------------------------------------------------------------------------------------------

Playmol is a molecular model builder which intends to facilitate the construction of initial setups
for Molecular Dynamics simulations of dense systems. It was initially designed as a bridge between
the molecular packing generator [Packmol] and the popular Molecular Dynamics simulator [LAMMPS].

The principle of Playmol consists in the bottom-up assembly of molecules and molecular systems by a
series of steps like:

1. __Typification__: atom types are defined and combinations of these types are registered as
templates for the identification of chemical bonds, angles, dihedrals, and impropers. These are the
topological features most usually considered in classical force fields, whose parameters can also be
fed to Playmol in this first step.

2. __Chemical topology definition__: atoms are created and explicitly connected via chemical bonds.
This is all Playmol needs to identify chemical compounds and automatically determine the angles and
proper dihedrals present in their molecules. Playmol does not automatically detect impropers, but
one can point them out manually.

3. __Instantiation__: chemical compounds identified in the previous step remain as abstract entities
until actual molecules are positioned in space. The [build] command makes it possible to instantiate
multiple molecules of each compound in order to create a molecular system.

4. __Replication and packing__: after instantiating at least one molecule of each compound, Playmol
can invoke [Packmol] to instantiate replicas in random positions, with possibly random orientations,
while trying to avoid undesirable overlaps. This allows a dense packing of molecules to be a viable
initial configuration for molecular dynamics.

5. __Formatting and storage__: finally, the system configuration and the force field parameters are
stored in a file whose format can be understood by a simulation software such as [LAMMPS].

Practical examples of the use of Playmol can be found in the section describing all Playmol
[commands] and in a special section with sample [scripts].


----------------------------------------------------------------------------------------------------
Playmol Scripts
----------------------------------------------------------------------------------------------------

Playmol interprets scripts composed of specific [commands]. Each command is a sequence of keywords
and parameter values separated by spaces and/or tabs. A script can include comments, identified by
the comment mark "#", which instructs Playmol to ignore all trailing characters in the same line. A
single command can span several lines by means of continuation marks "&", as in the example below
for the [bond] command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bond	C2       & # This is the central carbon atom,
    	C1  C3   & # which makes bonds with two other carbons
    	H21 H22    # and two hydrogens.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Except when a [for/next] construct or an [if/then/else] construct is found, the execution flow of a
Playmol script is sequential. Interpretation of every command line will always follow the three
steps below:

1. __Substitution of variables__: after reading a command line, Playmol will search for references
to [variables]. When a reference is found, it will check whether the variable has actually been
assigned by a previous [define] command. If so, Playmol will substitute the reference by the
corresponding value. This step is repeated until no references are found.

2. __Evaluation of mathematical expressions__: after replacing all variables, Playmol will look for
[mathematical/logical expressions](basics.html#math) in the command line, evaluate them, and substitute
every expression by its numerical result.

3. __Parsing and execution__: after finishing the two previous steps, Playmol will parse the command
line and try to execute it, in accordance with the syntaxes described in the [commands] section. The
program will post a message and stop if an error occurs.


----------------------------------------------------------------------------------------------------
<a name="variables"></a>
Variables
----------------------------------------------------------------------------------------------------

Variable assignment is done via the [define] command. For instance:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define compound as H2O
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following rules must be observed for choosing a variable name:

* It is case-sensitive.
* The first character must necessarily be a letter (A-Z or a-z).
* The remaining characters can be letters (A-Z or a-z), numbers (0-9), or underscores (_).

A variable can be referenced later on by either preceding its name with a symbol "$" or enclosing it
between symbols "${" and "}", with no spaces. For example:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write  xyz    $compound.xyz
write  lammps ${compound}_config.lmp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, Playmol would save an xyz file named `H2O.xyz` and a LAMMPS data file named
`H2O_config.lmp`. Note that the use of curly brackets is required in the second line because Playmol
would otherwise look for a variable named `compound_config`.


----------------------------------------------------------------------------------------------------
<a name="math"></a>
Mathematical/Logical Expressions
----------------------------------------------------------------------------------------------------

Playmol parses and evaluates mathematical or logical expressions placed inside curly brackets. This
is done once all variable references have been properly substituted. After evaluating, Playmol will
substitute the brackets and the expression by the obtained numerical result. For instance:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
define  factor  as 0.602214
define  density as { 1.3*$factor }
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mathematical expressions might result in an integer number or a real number. If it is a logical
expression, the result might be either an integer `1` (meaning `true`) or an integer `0` (meaning
`false`). Real numbers might be expressed using a decimal point (e.g. `2.0`) and/or scientific
notation (e.g. `1.38E-23`). The order of precedence of mathematical/logical operators is given in
[Table 1] and the supported mathematical functions are listed in [Table 2]. Note that there is no
specific operator for the logical `NOT` in [Table 1]. It is in fact implemented as the function
`not()` present in [Table 2].

<a name="Table_1"></a> **Table 1**: Precedence order of mathematical operators

| Operator                                 | Description                             | Precedence  |
|:----------------------------------------:|:----------------------------------------|:-----------:|
| \|                                       | Logical OR                              | 1 (lowest)  |
| \&                                       | Logical AND                             | 2           |
| < <br> <= <br> > <br> >= <br> == <br> != | Relational operators                    | 3           |
| + <br> -                                 | Addition<br>Subtraction                 | 4           |
| * <br> / <br> %                          | Multiplication<br>Division<br>Remainder | 5           |
| ^                                        | Exponentiation                          | 6           |
| `function( )`                            | Functions (see [Table 2])               | 7           |
| ( )                                      | Parentheses                             | 8 (highest) |

<a name="Table_2"></a> **Table 2**: Available mathematical functions

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
[variables]:		#variables

<!-- External links -->
[scripts]:		scripts.html
[commands]:		commands.html
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

