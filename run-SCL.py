#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, math
from SCL import *
from aiwlib.vec import *
#from aiwlib.iostream import *
#from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=(200., '#время релаксации'), t_eq=(10., '#время накопления равновесной статистики'), 
            steps=(10, '#число шагов между сбросом данных'), _repo='SCL', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            Nx=(32, '#размер решетки по X'), Ny=(32, '#размер решетки по y'), Nz=(32, '#размер решетки по Z'), spectr=(0, '#скважность сброса спектров (в больших шагах)'))
calc.tags.add('SC')
calc.mode = 'SCL'

model = calc.wrap(Model(), '', lambda: model.t/(calc.t_eq+calc.t_relax))
model.Hext[2] = calc.Hz
model._core.mesh_sz = ind(calc.Nx, calc.Ny, calc.Nz)
#print model.mesh_sz, model._core.mesh_sz

#if model.dense_mode==0 and (model.T>model.sparse*1.6 or model.T<model.sparse*1.6-.5): exit(0) 

model.init(calc.path)


#if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')
#fsph_av = File(calc.path+'f_av.sph', 'w')
istep = 0
while model.t<calc.t_relax+calc.t_eq:
    model.calc_eq = calc.t_eq and model.t>calc.t_relax
    model.calc(calc.steps)
    if calc.spectr and istep%calc.spectr==0:  model.calc_spectrum()
    istep += 1
    #if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()

model.finish()
