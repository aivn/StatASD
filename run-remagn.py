#!/usr/bin/python2
import os, math, shutil
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=50., t_remagn=100., steps=50, _repo='remagn', Hz=-0.5, M0z=1., mode=cvar.mode, f_rank=2, tau=0.)

model = calc.wrap(Model())
model.stoch0 = False
model.helic0 = False
model.entropy0 = False
model.helic_z = .1
model.parallel = True
model.f_rank = -1
model.data_rank = 5
model.fz_sz = 100
model.J = 1.
model.dt = 1e-2
model.T = 1.
model.K = 0.1
model.gamma = 1.
model.alpha = 0.1
#model.Hext = vecf(0., 0., calc.Hz)
model.nK = vecf(0., 0., 1.)
model.M0 = vecf(0., 0., calc.M0z)
model.patchT = False
model.T_sc_weight = .1

model.calc_eq = False

Tc = 1.5 if cvar.mode=='SC' else 2.2 if cvar.mode=='VCC' else 3.3
if Tc-.8<=model.T and model.T<=Tc+1.4: calc.t_relax *= 2
if Tc-.4<=model.T and model.T<=Tc+.7: calc.t_relax *= 2
if Tc-.2<=model.T and model.T<=Tc+.3: calc.t_relax *= 4

cond_name = 'init_conditions/%s-R%i-T%g-K%g.spins'%(cvar.mode, model.data_rank, model.T, model.K)
if os.path.exists(cond_name) and model.init(calc.path, cond_name): calc.t_max, t_init = calc.t_remagn, 0
else:
    model.init(calc.path)
    calc.t_max, t_init = calc.t_relax + calc.t_remagn, calc.t_relax
    print 'calc', cond_name
    while model.t<calc.t_relax: model.calc(calc.steps); calc.set_progress(model.t/calc.t_max, 'calc')
    model.dump_data(cond_name)
    shutil.move(calc.path+'tvals.dat', calc.path+'tvals0.dat')
    model.t = 0
    model.open_tvals(calc.path)
    
model.Hext[2] = calc.Hz

if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')

model.dump_Ms_arr()
while model.t<calc.t_remagn:
    model.calc(calc.steps); model.dump_Ms_arr()
    if not calc.tau and model.M[2]<0: calc.tau = model.t - calc.t_relax
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()
    calc.set_progress((t_init+model.t)/calc.t_max, 'calc')

model.finish()
