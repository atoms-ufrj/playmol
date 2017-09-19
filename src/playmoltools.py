#!/usr/bin/env python

import sys
import datetime
import argparse

#---------------------------------------------------------------------------------------------------

def amber2playmol( inp, out ):

  atoms = {}
  mass = {}
  note = {}
  params = {}

  bond = []
  angle = []
  dihedral = []
  improper = []
  fast = []
  equiv = []

  def prefix( types ):
    if all([word in mass for word in types.split()]):
      return ""
    else:
      return "#"    

  blocks = ['title','atom','hydrophilic','bond','angle','dihedral','improper',
            'fast','equivalence','type','params']
  m = 0
  for line in inp:
    string = line.rstrip()
    if (not string):
      m += 1
      continue
    elif (string == 'END'):
      break
    block = blocks[m]

    if (block == 'title'):
      title = string
      m += 1

    elif (block == 'hydrophilic'):
      hydrophilic = string
      m += 1

    elif (block == 'atom'):
      word = string.split()
      mass[word[0]] = word[1]
      note[word[0]] = ' '.join(word[3:])

    elif (block == 'bond'):
      types = string[0:5].replace('-',' ')
      bond.append(prefix(types)+'\t'.join(['bond_type',types,'$bond_style',string[5:22]]))

    elif (block == 'angle'):
      types = string[0:8].replace('-',' ')
      angle.append(prefix(types)+'\t'.join(['angle_type',types,'$angle_style',string[8:28]]))

    elif (block == 'dihedral'):
      types = string[0:11].replace('-',' ').replace('X','*')
      if not types.strip():
        types = dihedral[-1][14:25]
      IDIVF, PK, PHASE, PN = string[11:54].split()
      K = '{'+PK+'/'+IDIVF+'}'
      n = str(abs(int(float(PN))))
      d = str(int(float(PHASE)))
      dihedral.append(prefix(types)+'\t'.join(['dihedral_type',types,'$dihedral_style',K,n,d,'0']))

    elif (block == 'improper'):
      types = (string[6:9] + string[0:6] + string[9:11]).replace('-',' ').replace('X','*')
      PK, PHASE, PN = string[11:54].split()
      K = str(float(PK)/float(IDIVF))
      n = str(abs(int(float(PN))))
      if (abs(float(PHASE)) < 0.01):
        improper.append(prefix(types)+'\t'.join(['improper_type',types,'$improper_style',K,' 1',n]))
      elif (abs(float(PHASE) - 180.0) < 0.01):
        improper.append(prefix(types)+'\t'.join(['improper_type',types,'$improper_style',K,'-1',n]))
      else:
        exit('ERROR: unsupported improper type')

    elif (block == 'fast'):
      fast.append(string)

    elif (block == 'equiv'):
      equiv.append(string)

    elif (block == 'type'):
      typestr = string
      m += 1

    elif (block == 'params'):
      word = string.split()
      epsilon = word[2]
      sigma = str(round(1.78179743628*float(word[1]),5))
      params[word[0]] = '\t'.join([epsilon, sigma])

  print >>out, '#', title
  print >>out, '\n# Model definition:'
  print >>out, 'define\tpair_style     as lj/cut/coul/long'
  print >>out, 'define\tbond_style     as harmonic'
  print >>out, 'define\tangle_style    as harmonic'
  print >>out, 'define\tdihedral_style as charmm'
  print >>out, 'define\timproper_style as cvff'

  print >>out, '\n# Atom types:'
  for key, value in mass.iteritems():
    if key in params:
      print >>out, '\t'.join(['atom_type', key, '$pair_style', params[key], '# '+note[key]])
    else:
      print >>out, '\t'.join(['atom_type', key, '$pair_style', "UNDEFINED", '# '+note[key]])

  print >>out, '\n# Masses:'
  for key, value in mass.iteritems():
    print >>out, 'mass\t', key, '\t', value

  print >>out, '\n# Bond types:'
  for item in bond:
    print >>out, item

  print >>out, '\n# Angle types:'
  for item in angle:
    print >>out, item

  print >>out, '\n# Dihedral types:'
  for item in dihedral:
    print >>out, item

  print >>out, '\n# Improper types:'
  for item in improper:
    print >>out, item

#---------------------------------------------------------------------------------------------------

def pdb2playmol( inp, out, prepfile ):

  atoms = []
  bonds = []
  nameorder = []
  namecount = {}
  diameter = {}

  for line in inp:
    item = line.split()
    record = item[0]
    if (record == 'ATOM') or (record == 'HETATM'):

      atomType = item[-1]
      diameter[atomType] = str(2.0*float(item[-2]))

      atomName = item[2]
      if (atomName in namecount):
        namecount[atomName] += 1
      else:
        namecount[atomName] = 1

      item[1] = atomName + '_' + str(namecount[atomName])

      atoms.append(item[1:])

    elif (record == 'CONECT'):
      atom1 = item[1]
      i = int(atom1)-1
      for atom2 in item[2:]:
        j = int(atom2)-1
        if (([i,j] not in bonds) and ([j,i] not in bonds)):
          bonds.append([i,j])

  if (prepfile):
    prep = prep_dict( prepfile )

  print >>out, "# Atom definitions:"
  for atom in atoms:
    if namecount[atom[1]] > 1:
      atomname = atom[0]
    else:
      atomname = atom[1]
    if (prepfile):
      prepitem = prep[atom[2]][atom[1]]
      print >>out, '\t'.join(["atom", atomname, prepitem[0], prepitem[-1]])
    else:
      print >>out, '\t'.join(["atom", atomname, atom[-1]])

  print >>out, "\n# Bond definitions:"
  for bond in bonds:
    i, j = bond
    print >>out, '\t'.join(["bond", atoms[i][1], atoms[j][1]])

  print >>out, "\n# Molecular structure (xyz):"
  print >>out, "\nbuild"
  print >>out, len(atoms)
  for atom in atoms:
    print >>out, atom[1], atom[4], atom[5], atom[6]

#---------------------------------------------------------------------------------------------------

def pqr2playmol( inp, out ):

  atoms = []
  bonds = []
  nameorder = []
  namecount = {}
  diameter = {}

  for line in inp:
    item = line.split()
    record = item[0]
    if (record == 'ATOM') or (record == 'HETATM'):

      atomType = item[-1]
      diameter[atomType] = str(2.0*float(item[-2]))

      atomName = item[2]
      if (atomName in namecount):
        namecount[atomName] += 1
      else:
        namecount[atomName] = 1

      item[1] = atomName + '_' + str(namecount[atomName])

      atoms.append(item[1:])

    elif (record == 'CONECT'):
      atom1 = item[1]
      i = int(atom1)-1
      for atom2 in item[2:]:
        j = int(atom2)-1
        if (([i,j] not in bonds) and ([j,i] not in bonds)):
          bonds.append([i,j])

  print >>out, "# Atom definitions:"
  for atom in atoms:
    if namecount[atom[1]] > 1:
      atom[1] = atom[0]
    print >>out, '\t'.join(["atom", atom[1], atom[-1], atom[7]])

  print >>out, "\n# Bond definitions:"
  for bond in bonds:
    i, j = bond
    print >>out, '\t'.join(["bond", atoms[i][1], atoms[j][1]])

  print >>out, "\n# Molecular structure (xyz):"
  print >>out, "\nbuild"
  print >>out, len(atoms)
  for atom in atoms:
    print >>out, atom[1], atom[4], atom[5], atom[6]

  print >>out, "\n# Diameters for packing:"
  for key, value in diameter.iteritems():
    print >>out, '\t'.join(["diameter", key, value])


#---------------------------------------------------------------------------------------------------

def prep_dict( prepfile ):
  input = open(prepfile,'r')
  prep = {}
  waiting = True
  for line in input:
    word = line.split()
    if waiting and word:
      residue = word[0]
      prep[residue] = {}
      waiting = False
    elif word:
      if word[0] == "DONE":
        waiting = True
      elif word[0] == "STOP":
        return prep
      elif word[0].isdigit() and int(word[0]) > 3:
        prep[residue][word[1]] = word[2:]
  input.close()
  return prep

#---------------------------------------------------------------------------------------------------
# MAIN PROGRAM
#---------------------------------------------------------------------------------------------------

allowedFormats = ['amber','pdb','pqr']

formatDescriptions = """
Supported formats are:
   amber - AMBER force field parameter file
   pdb   - Protein Data Bank (PDB) file
   pqr   - Modified PDB file (with charges and radii)
"""

parser = argparse.ArgumentParser(
  description = 'Playmol Tools: convert files of different formats into Playmol scripts',
  formatter_class = argparse.RawTextHelpFormatter,
  epilog = formatDescriptions )

parser.add_argument('-f', metavar='format', dest='format', choices = allowedFormats,
                    required = True, help='input file format (required)')

parser.add_argument('-i', metavar='name',  dest='input',
                    help='input file name (default: stdin)')

parser.add_argument('-o', metavar='name', dest='output',
                    help='output file name (default: stdout)')

parser.add_argument('-p', metavar='name', dest='prepfile',
                    help='prep file name (default: None)')

args = parser.parse_args()

if (args.input):
  input = open(args.input,'r')
else:
  input = sys.stdin;

if (args.output):
  output = open(args.output,'w')
else:
  output = sys.stdout;

if (args.format == 'amber'):
  amber2playmol( input, output )
elif (args.format == 'pdb'):
  pdb2playmol( input, output, args.prepfile )
elif (args.format == 'pqr'):
  pqr2playmol( input, output )

print >>output, '\n# Generated by playmoltools on', datetime.datetime.today()

if (args.input):
  input.close()

if (args.output):
  output.close()
