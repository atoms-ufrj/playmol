/*
Language: Playmol
Description: a script-based software for building molecular models
Repository: https://github.com/atoms-ufrj/playmol
Author: Charlles Abreu <abreu@eq.ufrj.br>
Category: scientific
*/

function(hljs) {

  var VARIABLE = {
    variants: [
      { 
        begin: /\$\{[a-zA-Z]\w*\}/,
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
        begin: /\$[a-zA-Z]\w*/
      }
    ]
  };

  var NUMBER = {
    className: 'number',
    begin: /[-+]?[0-9]*\.?[0-9]+([dDeE][-+]?[0-9]+)?/
  };

  var KEYWORDS = {
    keyword:
      'define as for from in to downto next if then else endif atom_type mass ' +
      'diameter bond_type angle_type dihedral_type improper_type atom charge ' +
      'bond link build align include shell aspect packmol tolerance seed retry ' +
      'nloops fix copy pack diameter',
    built_in:
      'not abs exp log ln sqrt sinh cosh tanh sin cos tan asin acos atan ' +
      'int nint ceil floor mol'
  };

  var SPECIAL_KEYWORDS = {
    className: 'keyword',
    variants: [
      { begin: /\bbox\s+(lengths|density|volume|angles)\b/ },
      { begin: /\bwrite\s+(playmol|lammps|summary|xyz|lammpstrj)\b/ },
      { begin: /\baction\s+(execute|setup)\b/ },
      { begin: /\b(pre|suf)fix\s+(atom|type)s(\s+none)?\b/ },
      { begin: /\bimproper(\s+search)?\b/ },
      { begin: /\bextra\s+(bond|angle|dihedral)\b/ },
      { begin: /\bquit(\s+all)?\b/ },
      { begin: /\breset\s+all/ },
      { begin: /\breset(\s+((bond|angle|dihedral|improper)_types|atoms|charges|bonds|impropers|xyz|packmol))+/ }
    ]
  };

  return {
    aliases: ['playmol', 'pmol', 'mol'],
    keywords: KEYWORDS,
    contains: [
      hljs.inherit(hljs.APOS_STRING_MODE, {className: 'string', relevance: 0}),
      hljs.inherit(hljs.QUOTE_STRING_MODE, {className: 'string', relevance: 0}),
      hljs.COMMENT('#', '$', {relevance: 0}),
      hljs.C_NUMBER_MODE,
      VARIABLE,
      SPECIAL_KEYWORDS
    ]
  };
}
