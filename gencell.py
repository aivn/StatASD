#!/usr/bin/python3
'''usage: ./gencell.py coord1 coord2 coord2 ... ==> C++ code for crystal cell 
exmple for BCC: ./gencell.py 0,0,0 .5,.5,.5
'''

import os, sys
from aiwlib.vec import *

if len(sys.argv)==1: print(__doc__); exit(1)

atoms = [Vec(*map(float, l.split(','))) for l in sys.argv[1:]]
#-------------------------------------------------------------------------------
cell_sz = len(atoms); links, Qtable, Qtsz = [[] for i in range(cell_sz)], [[] for i in range(cell_sz)], []

for aID in range(cell_sz):
    l_min = 1e-8+min([(pos+atoms[aID2]-atoms[aID]).abs() for pos in Range(ind(-1,-1,-1), ind(2,2,2)) for aID2 in range(cell_sz) if not (pos==ind(0) and aID2==aID)])
    for pos in Range(ind(-1,-1,-1), ind(2,2,2)):
        for aID2 in range(cell_sz):
            if not (pos==ind(0) and aID2==aID) and (pos+atoms[aID2]-atoms[aID]).abs()<l_min: links[aID].append((ind(tuple(pos)), aID2))
    #print links[aID]
    for i, l_i in enumerate(links[aID]):
        for j, l_j in enumerate(links[aID][:i]): Qtable[aID].append(((l_i[0]+atoms[l_i[1]]-l_j[0]-atoms[l_j[1]]).abs(), i, j))
    Qtable[aID].sort(); Qtsz.append([]); l = 0 #Qtable[aID][0][0]
    for t in Qtable[aID]:
        if l!=t[0]: Qtsz[-1].append(1); l = t[0]
        else: Qtsz[-1][-1] += 1
#-------------------------------------------------------------------------------
nb_sz, Q_sz = len(links[0]), len(set(t[0] for t in Qtable[0]))
print('const int cell_sz = %i, nb_sz = %i, Q_sz = %i;'%(cell_sz, nb_sz, Q_sz))
print('const Link nb_pos[%i][%i] = { %s };'%(cell_sz, nb_sz,
                                             ',\n'.join('{ %s }'%', '.join('Link(%i, %i, %i, %i)'%(tuple(pos)+(k,)) for pos, k in L) for L in links)))
print('const Ind<2> Qtable[%i][%i] = { %s };'%(cell_sz, len(Qtable[0]),
                                             ',\n'.join('{ %s }'%', '.join('ind(%i, %i)'%t[1:] for t in L) for L in Qtable)))
print('const int Qtable_sz[%i][%i] = { %s };'%(cell_sz, Q_sz,
                                             ', '.join('{ %s }'%', '.join(map(str, L)) for L in Qtsz)))
#-------------------------------------------------------------------------------


          
