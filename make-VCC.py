#!/usr/bin/python2

from aiwlib.vec import *

atoms = [vec(0.,0.,0.), vec(.5,.5,.5)] # coords
# links by subcells
links = [[vec(0,0,0), vec(-1,0,0), vec(0,-1,0), vec(0,0,-1), vec(0,-1,-1), vec(-1,0,-1), vec(-1,-1,0), vec(-1,-1,-1)], [None]*8]

for i in range(8): links[1][i] = -links[0][i]
csphs = [[{} for j in range(8)] for sc in (0,1)]  # link length: k by sc and j

for sc in (0, 1):   #sublattices
    for j in range(8):
        cj = links[sc][j]  # cell
        rj = cj+atoms[1-sc] # coord
        for k in range(8):
            ck = links[1-sc][k]  
            rk = cj+ck+atoms[sc]
            d = int((rk-atoms[sc]).abs()*1000)
            if d: csphs[sc][j].setdefault(d, []).append(k)

for d in csphs[0]: print d
print '----'
for d in csphs[1]: print d
            
print 'const Link QW[2][8+3*8+3*8+8] = {'
for sc in (0, 1):   #sublattices
    print '    {'
    for j in range(8):
        print '    Link(%i, %i, %i, %i),'%(tuple(links[sc][j])+(1-sc,)),
        for d, Lk in sorted(csphs[sc][j].items()):
            for k in Lk: print 'Link(%i, %i, %i, %i),'%(tuple(links[1-sc][k])+(sc,)),
        print ''
    print '},'
print '};'
