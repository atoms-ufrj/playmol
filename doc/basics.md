Playmol Basics      {#basics}
==============

------------
Introduction
------------

Playmol scripts are text files containing commands described [here](@ref commands). Each command is a sequence of keywords and parameter values separated by spaces and/or tabs. A script can include comments, identified by the comment mark "#". When such mark is found, all trailing characters in the same line (including the comment mark itself) are ignored. A single command can span several lines by means of continuation marks "&". When the actual part of a command line (i.e., excluding comments) ends with a continuation mark, the command will continue in the next non-empty line.

In all examples described in the [commands](@ref commands) section, the units employed for physically meaningful values are those corresponding to [LAMMPS real units].

---------
Variables
---------

------------------------
Mathematical Expressions
------------------------


