#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from LL3 import *
from aiwlib.vec import *
from aiwlib.iostream import *
from aiwlib.SphereF import *
from aiwlib.racs import *

calc = Calc(t_relax=(50., '#время релаксации'), t_eq=(10., '#время накопления равновесной статистики'),
            steps=(10, '#число шагов между сбросом данных'), _repo='chi', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            dH=(.01, '#приращение внешнего поля для расчета восприимчивости'),
            Tmax=(6., '#максимальная температура'), Tmin=(.1, '#минимальная температура'), dT=(.1, '#шаг температуры'), mode=cvar.mode,
            init=('', '#путь к прерванному ранее расчету'))
calc.tags.add(cvar.mode)

Tc = 1.445 if cvar.mode=='SC' else 2.05 if cvar.mode=='BCC' else 3.17

model = calc.wrap(Model(), '', lambda: (calc.Tmax-model.T)/(calc.Tmax-calc.Tmin))
model.M0[2] = 0
model.Hext[2] = calc.Hz
model._core.T = calc.Tmax
model.corr_max = 0

model.init(calc.path)
if calc.init:
    model.load_data(calc.init+'/spins')
    os.system('cp %s/chi.dat %s/chi.dat'%(calc.init, calc.path))
    fchi = open(calc.path+'chi.dat', 'a')
    #racs = pickle.load(open(calc.init+'/.RACS'))
    model._core.T = float(open('%s/chi.dat'%calc.init).readlines()[-1].split()[0])
    fcontinue = True
else:
    fchi = open(calc.path+'chi.dat', 'w')
    fchi.write('#:T chi M0 M1 W0 W1 eta0 eta1\n')
    fcontinue = False

while model.T>=calc.Tmin:
    t_relax, dT = calc.t_relax, calc.dT
    if Tc-.8<=model.T and model.T<=Tc+1.4: t_relax *= 2; #dT /= 2
    if Tc-.4<=model.T and model.T<=Tc+.7:  t_relax *= 4; dT /= 2
    if Tc-.2<=model.T and model.T<=Tc+.3:  t_relax *= 8; dT /= 2
    if fcontinue: model.T -= dT; fcontinue = False

    t_max0 = model.t + t_relax; t_max1 = t_max0 + calc.t_eq
    model.clean_av_eq(); model.Hext[2] = calc.Hz+calc.dH
    while model.t<t_max1: model.calc_eq = model.t>t_max0; model.calc(calc.steps)
    model.finish(); M1, W1, eta1 = model.Meq.abs(), model.Weq[0], model.eta_eq[0]

    t_max0 = model.t + t_relax; t_max1 = t_max0 + calc.t_eq
    model.clean_av_eq(); model.Hext[2] = calc.Hz
    while model.t<t_max1: model.calc_eq = model.t>t_max0; model.calc(calc.steps)
    model.finish(); M0, W0, eta0 = model.Meq.abs(), model.Weq[0], model.eta_eq[0]

    print>>fchi, model.T, (M1-M0)/calc.dH, M0, M1, W0, W1, eta0, eta1; fchi.flush()
    model.dump_data(calc.path+'spins')
    model._core.T -= dT
    
    
