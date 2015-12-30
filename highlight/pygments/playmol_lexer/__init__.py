# -*- coding: utf-8 -*-
"""
    pygments.lexers.playmol
    ~~~~~~~~~~~~~~~~~~~~~~~

    Lexer for Playmol.

    :copyright: Copyright 2006-2013 by the Pygments team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from pygments.lexer import Lexer, RegexLexer, include, bygroups, using, \
     this, combined, inherit, do_insertions
from pygments.token import Text, Comment, Operator, Keyword, Name, String, \
     Number, Punctuation, Error, Literal, Generic
from pygments.scanner import Scanner

__all__ = ['PlaymolLexer']

class PlaymolLexer(RegexLexer):
    """
    Lexer for Playmol.

    """
    name = 'Playmol'
    aliases = ['playmol']
    filenames = ['*.mol', '*.pmol', '*.playmol']
    mimetypes = ['text/plain']

    tokens = {
        'root': [
            (r'#.*\n', Comment),
            include('strings'),
            include('core'),
            (r'\$[a-z][a-z0-9_]*', Name.Variable),
            (r'\$\{[a-z][a-z0-9_]*\}', Name.Variable),
            include('nums'),
        ],
        'core': [
            # Statements
            (r'\b(atom_type|mass|bond_type|angle_type|dihedral_type|improper_type|atom|'
             r'charge|bond|link|extra\s+(bond|angle|dihedral)|improper\s+search|improper|'
             r'xyz|build|align|box\s+(lengths|density|volume|angles)|packmol|'
             r'write\s+(playmol|lammps|summary|xyz|lammpstrj)|prefix\s+(atoms|types)|'
             r'prefix|suffix\s+(atoms|types)|suffix|to|downto|next|if|then|else|endif|'
             r'include|reset|clean_types|shell|quit all|quit|tolerance|seed|retry|nloops|'
             r'fix|copy|pack|action\s+(execute|setup)|aspect|)\s*\b',
             Keyword),

            # Operators
            (r'(\^|\*|\+|-|\/|\%|<|>|<=|>=|==|!=)', Operator),

#            (r'\{()\}', Punctuation),

            # Intrinsics
            (r'\b(not|abs|exp|log|ln|sqrt|sinh|cosh|tanh|sin|cos|'
             r'tan|asin|acos|atan|int|nint|ceil|floor)\s*\b',
             Name.Builtin),

            # Comparing Operators
            (r'(&|\|)', Operator.Word),
        ],

        'strings': [
            (r'(?s)"(\\\\|\\[0-7]+|\\.|[^"\\])*"', String.Double),
            (r"(?s)'(\\\\|\\[0-7]+|\\.|[^'\\])*'", String.Single),
        ],

        'nums': [
            (r'\d+(?![.Ee])', Number.Integer),
            (r'[+-]?\d*\.\d+([eE][-+]?\d+)?', Number.Float),
            (r'[+-]?\d+\.\d*([eE][-+]?\d+)?', Number.Float),
        ],
    }

