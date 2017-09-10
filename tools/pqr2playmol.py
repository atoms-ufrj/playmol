#!/usr/bin/env python

import fileinput

atoms = []
bonds = []
nameorder = []
namecount = {}
diameter = {}

for line in fileinput.input():
    item = line.split()
    record = item[0]
    if (record == 'ATOM') or (record == 'HETATM'):

        atomType = item[-1]
        diameter[atomType] = item[-2]

        atomName = item[2]
        if (atomName in namecount):
            namecount[atomName] += 1
        else:
            namecount[atomName] = 1

        if (namecount[atomName] > 1):
            item[1] = atomName + str(namecount[atomName])
        else:
            item[1] = atomName

        atoms.append(item[1:])
    elif (record == 'CONECT'):
        atom1 = item[1]
        i = int(atom1)-1
        for atom2 in item[2:]:
            j = int(atom2)-1
            if (([i,j] not in bonds) and ([j,i] not in bonds)):
                bonds.append([i,j])

print "# Atom definitions:"
for atom in atoms:
    print '\t'.join(["atom", atom[0], atom[-1], atom[7]])

print "\n# Bond definitions:"
for bond in bonds:
    i, j = bond
    print '\t'.join(["bond", atoms[i][0], atoms[j][0]])

print "\n# Molecular structure (xyz):"
print "\nbuild"
print len(atoms)
for atom in atoms:
    print atom[0], atom[4], atom[5], atom[6]

print "\n# Diameters for packing:"
for key, value in diameter.iteritems():
    print '\t'.join(["diameter", key, value])
