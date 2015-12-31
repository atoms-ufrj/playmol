Playmol Basics      {#basics}
==============

Introduction
------------

Playmol scripts are text files containing commands described [here](commands.html). Each command is a sequence of keywords and parameter values separated by spaces and/or tabs. A script can include comments, identified by the comment mark "#". When such mark is found, all trailing characters in the same line (including the comment mark itself) are ignored. A single command can span several lines by means of continuation marks "&". When the actual part of a command line (i.e., excluding comments) ends with a continuation mark, the command will continue in the next non-empty line.

In all examples described in the [commands] section, the units employed for physically meaningful values are those corresponding to [LAMMPS real units].

---------
Variables
---------

------------------------
Mathematical Expressions
------------------------

Operator precedence table:

| Operator                                 | Description                               | Prececence  |
|:----------------------------------------:|:------------------------------------------|:-----------:|
| \|                                       | Logical OR                                | 1 - lowest  |
| \&                                       | Logical AND                               | 2           |
| < <br> <= <br> > <br> >= <br> == <br> != | Relational operators                      | 3           |
| + <br> -                                 | Addition<br> Subtraction                  | 4           |
| * <br> / <br> %                          | Multiplication<br> Division<br> Remainder | 5           |
| ^                                        | Exponentiation                            | 6           |
| _function_( )                            | Mathematical functions                    | 7           |
| ( )                                      | Parentheses                               | 8 - highest |

__Note__: There is not a specific operator for logical NOT. Its effect is obtained by using the function `not()`.

Available mathematical functions:

| Function | Result type      | Description                |
|:--------:|:----------------:|:---------------------------|
| abs      | Same as argument | Absolute value             |
| acos     | Real             | Inverse cosine             |
| asin     | Real             | Inverse sine               |
| atan     | Real             | Inverse tangent            |
| ceil     | Integer          | Largest previous integer   |
| cos      | Real             | Cosine                     |
| cosh     | Real             | Hyperbolic cosine          |
| exp      | Real             | Exponential                |
| floor    | Integer          | Smallest following integer |
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





[LAMMPS real units]:	http://lammps.sandia.gov/doc/units.html
[read_data]:		http://lammps.sandia.gov/doc/read_data.html
[xyz file format]:	http://openbabel.org/wiki/XYZ_(format)
[Packmol User's Guide]:	http://www.ime.unicamp.br/~martinez/packmol/quickguide/
[VMD]:			http://www.ks.uiuc.edu/Research/vmd/
[commands]:             commands.html

