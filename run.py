#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=(50., '#начальное время релаксации'), t_eq=(10., '#время накопления равновесной статистики'), t_max=(0., '#максимальное время расчета'),
            steps=(10, '#число шагов между сбросом данных'), _repo='repo', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            M0z=(1., '#начальная намагниченность (для задания из командной строки)'), spectr=(0, '# скважность сброса спектра (в больших шагах)'),
            _queue_main_par='mode', mode=cvar.mode, Tc=(0., '#оценка Tc для автоподбора t_relax'), dump_spins=(0, '#скважность сброса spins для продолжения расчета'))
calc.tags.add(cvar.mode)

model = calc.wrap(Model(), '', lambda: model.t/calc.t_max)

#model.M0[2] = calc.M0z
#model.Hext[2] = calc.Hz

#model.cL = dict((int(l.split()[0]), float(l.split()[1])) for l in open('LcL-v0.dat') if l.strip() and l[0]!='#')[1<<model.data_rank]

model.M0[2] = calc.M0z
model.init(calc.path)

if os.path.exists(calc.path+'spins'):
    model._core.t = float(open(calc.path+'tvals.dat').readlines()[-1].split()[0])
    if not model.load_data(calc.path+'spins'): print 'LOAD SPINS FAILED'; exit()
else:    
    model.Hext[2] = calc.Hz

    if calc.Tc: Tc = calc.Tc
    else: Tc = 1.45 if cvar.mode=='SC' else 2.1 if cvar.mode=='BCC' else 3.2

    calc.t_relax0 = calc.t_relax
    if not calc.t_max:
        if Tc-.8<=model.T and model.T<=Tc+1.4: calc.t_relax *= 2
        if Tc-.4<=model.T and model.T<=Tc+.7: calc.t_relax *= 2
        if Tc-.2<=model.T and model.T<=Tc+.3: calc.t_relax *= 4
        calc.t_max = calc.t_relax+calc.t_eq


if model.f_rank>=0: fsph = File(calc.path+'f.sph', 'a')
#fsph_av = File(calc.path+'f_av.sph', 'w')

istep = 0

while model.t<calc.t_max:
    model.calc_eq = calc.t_eq and model.t>calc.t_relax
    model.calc(calc.steps)    
    if model.f_rank>=0:  model.f.dump(fsph); fsph.flush()
    if calc.spectr and istep%calc.spectr==0: model.calc_spectrum()
    if calc.dump_spins and istep%calc.dump_spins==0: model.dump_data(calc.path+'spins')
    istep += 1
    
model.finish()
