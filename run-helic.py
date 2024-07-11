#!/usr/bin/python2
import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(steps=10, _repo='repo', Hz=0., M0z=1., mode=cvar.mode)

model = calc.wrap(Model())
model.stoch0 = False
model.helic0 = not model.stoch0
model.entropy0 = False
model.helic_z = .1
model.helic_n = 1
model.parallel = True
model.f_rank = 2
model.f_use = False
model.data_rank = 5
model.J = 1.
model.dt = 1e-2
model.T = 1.
model.K = 0.
model.gamma = 1.
model.alpha = 0.1
model.Hext = vecf(0., 0., calc.Hz)
model.nK = vecf(0., 0., 1.)
model.M0 = vecf(0., 0., calc.M0z)
#model.mg_split = 1
#model.T = .0015*2**model.mg_split
#model.zeta = 0.
model.cL = dict((int(l.split()[0]), float(l.split()[1])) for l in open('LcL-v0.dat') if l.strip() and l[0]!='#')[1<<model.data_rank]

try: Meq = dict(map(float, l.split()) for l in open('dat/SC-TM.dat') if l[0]!='#' and l.strip())[model.T]
except: print 'T not in dat/SC-TM.dat'; exit(1)

model.calc_eq = False
model.init(calc.path)

if model.f_use: fsph = File(calc.path+'f.sph', 'w')
#fsph_av = File(calc.path+'f_av.sph', 'w')

while model.M.abs()<Meq: 
    model.calc(calc.steps)
    
    if model.f_use: model.calc_f(); model.f.dump(fsph); fsph.flush()

    calc.set_progress(model.M.abs()/Meq, 'calc')

model.finish()
