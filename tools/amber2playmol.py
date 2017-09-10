#!/usr/bin/env python

import fileinput
import datetime

atoms = {}
mass = {}
params = {}

bond = []
angle = []
dihedral = []
improper = []
fast = []
equiv = []

blocks = ['title','atom','hydrophilic','bond','angle','dihedral','improper','fast','equivalence','type','params']
m = 0;
for line in fileinput.input():
    string = line.strip()
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

    elif (block == 'bond'):
        types = string[0:5].replace('-',' ')
        bond.append(' '.join(['bond_type\t',types,'\t$bond_style',string[5:22]]))

    elif (block == 'angle'):
        types = string[0:8].replace('-',' ')
        angle.append(' '.join(['angle_type\t',types,'\t$angle_style',string[8:28]]))

    elif (block == 'dihedral'):
        types = string[0:11].replace('-',' ').replace('X','*')
        IDIVF, PK, PHASE, PN = string[11:54].split()
        K = str(float(PK)/float(IDIVF))
        n = str(abs(int(float(PN))))
        d = str(float(PHASE))
        dihedral.append(' '.join(['dihedral_type\t',types,'\t$dihedral_style',K,n,d,'1.0']))

    elif (block == 'improper'):
        types = string[0:11].replace('-',' ').replace('X','*')
        PK, PHASE, PN = string[11:54].split()
        K = str(float(PK)/float(IDIVF))
        n = str(abs(int(float(PN))))
        if (float(PHASE) == 0.0):
            improper.append(' '.join(['improper_type\t',types,'\t$improper_style',K,'1',n]))
        elif (float(PHASE) == 180.0):
            improper.append(' '.join(['improper_type\t',types,'\t$improper_style',K,'-1',n]))
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
        params[word[0]] = word[1] + ' ' + word[2]

print '#', title
print
print '# File generated on ', datetime.date.today()
print
print '# Model definition:'
print 'define\tpair_style     as lj/cut/coul/long'
print 'define\tbond_style     as harmonic'
print 'define\tangle_style    as harmonic'
print 'define\tdihedral_style as charmm'
print 'define\timproper_style as cvff'

print '\n# Atom types:'
for key, value in params.iteritems():
    print 'atom_type\t', key, '\t$pair_style', value

print '\n# Masses:'
for key, value in mass.iteritems():
    print 'mass\t', key, '\t', value

print '\n# Bond types:'
for item in bond:
    print item

print '\n# Angle types:'
for item in angle:
    print item

print '\n# Dihedral types:'
for item in dihedral:
    print item

print '\n# Improper types:'
for item in improper:
    print item

