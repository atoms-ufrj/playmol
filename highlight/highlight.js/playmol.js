/*
Language: Playmol
Author: Charlles Abreu <abreu@eq.ufrj.br>
Category: scientific
*/

function(hljs) {

  var VARIABLE = {
    variants: [
      { 
        begin: /\$\{[a-zA-Z][a-zA-Z0-9_]*\}/,
        returnBegin: true,
        contains: [
          {
            className: 'variable',
            begin: /\$\{/,
            end: /\}/,
            excludeBegin: true,
            excludeEnd: true
          }
        ]
      },
      {
        className: 'variable',
        begin: /\$[a-zA-Z][a-zA-Z0-9_]*/
      }
    ]
  };

  var NUMBER = {
    className: 'number',
    begin: '(?=\\b|\\+|\\-|\\.)(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*)(?:[de][+-]?\\d+)?\\b\\.?'
  };

  var KEYWORDS = {
    keyword:
      'define as for from in to downto next if then else endif atom_type mass ' +
      'bond_type angle_type dihedral_type improper_type atom charge bond link ' +
      'build align include reset clean_types shell aspect ' +
      'packmol tolerance seed retry nloops fix copy pack',
    built_in:
      'not abs exp log ln sqrt sinh cosh tanh sin cos tan asin acos atan int nint ceil floor'
  };

  var EXTRA_KEYWORDS = {
    className: 'keyword',
    variants: [
      { begin: /^\s*box\s+(lengths|density|volume|angles)\b/ },
      { begin: /^\s*write\s+(playmol|lammps|summary|xyz|lammpstrj)\b/ },
      { begin: /\s+action\s+(execute|setup)\b/ },
      { begin: /^\s*(pre|suf)fix\s+(atom|type)s(\s+none)?\b/ },
      { begin: /^\s*improper(\s+search)?\b/ },
      { begin: /^\s*extra\s+(bond|angle|dihedral)\b/ },
      { begin: /^\s*quit(\s+all)?\b/ }
    ]
  };

  return {
    aliases: ['playmol', 'pmol', 'mol'],
    keywords: KEYWORDS,
    contains: [
      hljs.inherit(hljs.APOS_STRING_MODE, {className: 'string', relevance: 0}),
      hljs.inherit(hljs.QUOTE_STRING_MODE, {className: 'string', relevance: 0}),
      hljs.COMMENT('#', '$', {relevance: 0}),
      NUMBER,
      VARIABLE,
      EXTRA_KEYWORDS
    ]
  };
}
