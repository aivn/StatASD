#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=(30., '#время релаксации'), T0=(.1, '#начальная температура'), T1=(1.5, '#конечная температура'), dT=(.1, '#шаг температуры'),
            steps=(10, '#число шагов между сбросом данных'), _repo='repo', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            M0z=(1., '#начальная намагниченность (для задания из командной строки)'), mode=cvar.mode)
calc.tags.add(cvar.mode)

model = calc.wrap(Model(), '', lambda: model.t/t_max_tot)
model.M0[2] = calc.M0z
model.Hext[2] = calc.Hz

Tc = 1.5 if cvar.mode=='SC' else 2.2 if cvar.mode=='BCC' else 3.3

Ttr, t_max_tot = [], 0
T = calc.T0
while T<calc.T1:
    t_relax = calc.t_relax
    if Tc-.8<=T and T<=Tc+1.4: t_relax *= 2
    if Tc-.4<=T and T<=Tc+.7: t_relax *= 2
    if Tc-.2<=T and T<=Tc+.3: t_relax *= 2
    Ttr.append((T, t_relax))
    T += calc.dT
    t_max_tot += t_relax
    
model.init(calc.path)
print 'init OK'

for T, t_relax in Ttr:
    model._core.T, t_max = T, model.t + t_relax
    while model.t<t_max: model.calc(calc.steps)    
    model.dump_data('init_conditions/%s-R%i-T%g-H%g-K%g.spins'%(cvar.mode, model.data_rank, model.T, calc.Hz, model.K))

