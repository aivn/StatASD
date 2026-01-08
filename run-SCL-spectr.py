#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from SCL import *
from aiwlib.vec import *
#from aiwlib.iostream import *
#from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=(200., '#время релаксации'), t_eq=(0., '#время накопления равновесной статистики'), 
            steps=(100, '#число шагов между сбросом данных'), _repo='SCL-spectr', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            Nx=(32, '#размер решетки по X'), Ny=(32, '#размер решетки по y'), Nz=(32, '#размер решетки по Z'))
calc.tags.add('SC')
calc.mode = 'SCL'

model = calc.wrap(Model(), '', lambda: model.t/(2*calc.t_relax))
model.Hext[2] = calc.Hz
model._core.mesh_sz = ind(calc.Nx, calc.Ny, calc.Nz)
#print model.mesh_sz, model._core.mesh_sz
calc.T0 = model.T

#if model.dense_mode==0 and (model.T>model.sparse*1.6 or model.T<model.sparse*1.6-.5): exit(0) 

model.init(calc.path)
model.calc_eq = False


#if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')
#fsph_av = File(calc.path+'f_av.sph', 'w')

while model.t<2*calc.t_relax:
    if model.t>calc.t_relax:
        if model.T: model.calc_spectrum(calc.path+'spectrum1.dat')    
        model._core.T = 0 
    model.calc(calc.steps)
    model.calc_spectrum()
    #if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()

model._core.T = calc.T0
model.calc_spectrum(calc.path+'spectrum2.dat')    
#model.finish()
