#!/usr/bin/python2 
import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=50., t_eq=10., steps=10, _repo='repo', Hz=0., M0z=1., mode=cvar.mode)

model = calc.wrap(Model())
model.stoch0 = False
model.helic0 = False
model.entropy0 = False
model.helic_z = .1
model.parallel = True
model.f_rank = 2
model.data_rank = 5
model.fz_sz = 100
model.J = 1.
model.dt = 1e-2
model.T = 1.
model.K = 0.
model.gamma = 1.
model.alpha = 0.1
model.Hext = vecf(0., 0., calc.Hz)
model.nK = vecf(0., 0., 1.)
model.M0 = vecf(0., 0., calc.M0z)
model.patchT = False
model.T_sc_weight = 2.

#model.mg_split = 1
#model.T = .0015*2**model.mg_split
#model.zeta = 0.
#model.cL = dict((int(l.split()[0]), float(l.split()[1])) for l in open('LcL-v0.dat') if l.strip() and l[0]!='#')[1<<model.data_rank]

Tc = 1.5 if cvar.mode=='SC' else 2.2 if cvar.mode=='VCC' else 3.3
if Tc-.8<=model.T and model.T<=Tc+1.4: calc.t_relax *= 2
if Tc-.4<=model.T and model.T<=Tc+.7: calc.t_relax *= 2
if Tc-.2<=model.T and model.T<=Tc+.3: calc.t_relax *= 4
calc.t_max = calc.t_relax+calc.t_eq

model.init(calc.path)

if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')
#fsph_av = File(calc.path+'f_av.sph', 'w')

while model.t<calc.t_max:
    model.calc_eq = model.t>calc.t_relax
    model.calc(calc.steps)
    
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()

    calc.set_progress(model.t/calc.t_max, 'calc')

model.finish()
