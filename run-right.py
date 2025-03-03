#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_max=(100., '#максимальное время расчета'), t_right=(0., '#интервал между расчетами правой части, если !=0 заменяет dM'),
            steps=(10, '#число шагов между сбросом данных'), _repo='right', Hz=(-.1, '#внешнее поле (для задания из командной строки)'),
            dM=(.02, '#шаг намагниченности с которым считаются правые части'), tau=(0, '#время перемагничивания'), mode=cvar.mode)
calc.tags.add(cvar.mode)
calc.tags.add('right')

M0 = 1
model = calc.wrap(Model(), '', lambda: max(model.t/calc.t_max, (M0-model._core.M[2])/(2*M0)))
model.Hext[2] = 0
model.data_rank = 8
model.corr_max = 0


init_cond = 'init_conditions/%s-R%i-T%g-H%g-K%g.spins'%(cvar.mode, model.data_rank, model.T, 0, model.K)
if not os.path.exists(init_cond): print init_cond, 'not found'; exit(1)
model.Hext[2] = calc.Hz

model.init(calc.path)
if not model.load_data(init_cond): print init_cond, 'not loaded'; exit(1)
M0 = M1 = model.M[2]; t1 = 0
H, T, K = model.Hext[2], model.T, model.K
n_b = {'SC':6, 'BCC':8, 'FCC':12}[cvar.mode]


if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'w')
#fsph_av = File(calc.path+'f_av.sph', 'w')

fright = open(calc.path+'right.dat', 'w'); fright.write('#:t M eta UpsilonM Q\n')
gadt2 = 2*model.alpha*model.gamma*model.dt

while model.t<calc.t_max and model.M[2]>-M0:    
    model.calc(calc.steps)    
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()
    if model.M[2]<0 and not calc.tau: calc.tau = model.t
    if (calc.t_right and model.t>=t1+calc.t_right) or (not calc.t_right and model.M[2]<M1-calc.dM):
        model._core.T, model._core.Hext[2], model._core.K = 0., 0., 0.
        Mz1, eta1 = model.M[2], model.eta[0]
        model.calc(1)
        Mz2, eta2 = model.M[2], model.eta[0]
        model._core.T, model._core.Hext[2], model._core.K = T, H, K
        print>>fright, model.t-model.dt/2, (Mz1+Mz2)/2, (eta1+eta2)/2, (Mz2-Mz1)/(gadt2*n_b), (eta2-eta1)/(2*gadt2); fright.flush()
        if calc.t_right: t1 = model.t
        else: M1 = Mz2
        
