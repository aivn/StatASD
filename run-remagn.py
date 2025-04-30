#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math, shutil
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_max=(500., '# максимальное время перемагничивания'), steps=10, _repo='remagn', Hz=-0.1, mode=cvar.mode, tau_r=(0., '# время перемагничивания'))

M0 = 1
model = calc.wrap(Model(), '', lambda: max(model.t/calc.t_max, 1-model.M[2]/M0))
model.data_rank = 6
model.f_rank = -1
model.fz_sz = 0
model.calc_eq = False

model.init(calc.path)

cond_name = 'init_conditions/%s-R%i-T%g-H0-K%g.spins'%(cvar.mode, model.data_rank, model.T, model.K)
if not os.path.exists(cond_name) or not model.load_data(cond_name): print 'init cond. %r not found'%cond_name; exit(1)
print 1

calc.M0 = M0 = model.M[2]

model.Hext[2] = calc.Hz

if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')

#model.dump_Ms_arr()
while model.t<calc.t_max:
    model.calc(calc.steps); #model.dump_Ms_arr()
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()
    if model.M[2]<0: calc.tau_r = model.t; break

#model.finish()
