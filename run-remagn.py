#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_max=(500., '# максимальное время перемагничивания'), steps=10, _repo='remagn', Hz=(-0.1, '#амплитуда поля перемагничивания'),
            mode=cvar.mode, tau_r=(0., '# время перемагничивания'),
            angle=(0., '#угол между nK и Hext в градусах'))
calc.tags.add(cvar.mode)

M0 = 1
model = calc.wrap(Model(), '', lambda: max(model.t/calc.t_max, 1-model.M[2]/M0))
model.data_rank = 6
model.f_rank = -1
model.fz_sz = 0
model.calc_eq = False

model.init(calc.path)

cond_name = 'init_conditions/%s-R%i-T%g-H0-K%g.spins'%(cvar.mode, model.data_rank, model.T, model.K)
if not os.path.exists(cond_name) or not model.load_data(cond_name): print 'init cond. %r not found'%cond_name; exit(1)

calc.Mbeg = M0 = model.M[2]

model.Hext[0], model.Hext[2] = calc.Hz*math.sin(calc.angle*math.pi/180), calc.Hz*math.cos(calc.angle*math.pi/180)

if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')

#model.dump_Ms_arr()
while model.t<calc.t_max:
    model.calc(calc.steps); #model.dump_Ms_arr()
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()
    if model.fz_sz>=0:  model.dump_fz(calc.path+'/fz-t%06g.dat'%model.t, False)
    if model.f_eta_sz>=0:  model.dump_f_eta(calc.path+'/feta-t%06g.dat'%model.t, False)
    if model.M[2]<0: calc.tau_r = model.t; break

model.finish()

calc.Mend = model.M[2]
if model.fz_sz>=0:  model.dump_fz(calc.path+'/fz-eq.dat', True)
if model.f_eta_sz>=0:  model.dump_f_eta(calc.path+'/feta-eq.dat', True)

