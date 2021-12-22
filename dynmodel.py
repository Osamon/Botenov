# -*- coding: utf-8 -*-

"""
Фит для динамической модели
"""

# -*- coding: utf-8 -*-

import os
import argparse
import itertools

import h5py as h5
import numpy as np
from sklearn.linear_model import Lasso

from mc2py.scutil import sv
from mc2py.scutil import ld
from mc2py.h52obj import dump
from mc2py.zgiw import giwstate, SetMFADrive, step, getv, IcamRD, BurnTo, GetHgrp3, SetHgrp3, RelaxCore, set_upz

from dynkoeffs import get_koeffs
from bgk_koeffs import bgk_koeffs_XeW,translated_bgk
from initpoint_calcdims import set_initial_state

lambda_id = 2.895e-05
lambda_xe = 2.12e-05


def set_state_for_dyn(w_cur, NH, Hnums):
    v = getv()

    v.YMFLAGSTAT = 1
    v.step()
    v.YZBORMODE = 1
    v.step()
    #v.YZBORMODE = 2
    #v.step()

    SetHgrp3([1, 1, 1], NH, Hnums)
    for i in range(20):
        v.step()

    v.YMINTPOW_SET = w_cur
    v.YM_STROBE_POWSET = 1
    v.step(n=5)

    v.YMFLAGSTAT = 0
    v.step()
    v.YZBORMODE = 2
    v.YMINTPOW_SET_FLAG = 1
    for i in range(20):
        v.step()

    v.YZBORMODE = 1
    v.step()
    v.YZBORMODE = 2
    v.step()
    v["#YM#YMTIME"] = 0
    v.step()


def make_dyn_model(bgk, hdf_in):
    #Выбор данных из уже набранных файлов bgk_{}.h5
    print(hdf_in)
    for hdf in hdf_in:
        print(hdf)
        with h5.File(hdf, 'r') as f5:
            #f5=translated_bgk(f5)          
            LIcam = np.array(f5['I'])
            #LIcam = np.array([i for i in f5['I']])
            LHgrp = np.array(f5['h'])
            #LHgrp = np.array([i for i in f5['h']])
            wt = np.array(f5['wt'])
            offset = np.array(f5['offset'])
            LXe=np.array([bgk_koeffs_XeW(bgk,['Xe'],i)['Xe'] for i in f5['Xe']])
            LW=np.array([bgk_koeffs_XeW(bgk,['w'],i)['w'] for i in f5['w']])
            
    res = {'LIcam': LIcam, 'LHgrp': LHgrp, 'LXe': LXe, 'LW': LW, 'wt': wt, 'offset': offset}
    return res



def make_dyn_model_old(bgk, dims, defaults, initpoint, burn, norm_koef, plant, block, kamp):
    #Набор данных по модели
    v = getv()
    v.syncro(1)

    v["#YM#YMTIME"] = 0

    SetMFADrive(plant=plant, block=block, kamp=kamp)

    NH = dims['NH']
    Hnums = dims['Hnums']
    erheat = defaults['erheat'][-1]
    wnom = defaults['Wnom'][-1]

    set_initial_state(NH)

    #v.YQ_AKNP_TYPE = 1
    v.step()

    Icam0 = initpoint['I_nom']

    w_cur_max = 105 if burn < 320 else 100

    gs_init = giwstate()
    gs = giwstate()

    gs_init.fix()

    if burn:
        BurnTo(burn)
    gs.fix()

    v["#YM#YMTIME"] = 0

    xe_norm = norm_koef['xe_norm'][-1]

    SetMFADrive(W=100, hgrp=np.ones(NH)*1.03)
    v.YMFLAGSTAT = 1
    v.YZBORMODE = 2
    v.YMFLAG_SMCALC = 1
    #v.YQ_AKNP_TYPE = 1
    v.step()

    print('Burn {}'.format(burn))
    transient = []

    
    for w_cur in np.linspace(10, w_cur_max, 6):
        print('W {}'.format(w_cur))
        set_state_for_dyn(w_cur, NH, Hnums)

        if float(v.YMBOR_COR) > 0.5:
            hgrp = np.ones(len(Hnums))*1.03
            set_upz()
            if plant == 'kaln':
                hgrp[0] = 0.5 #В Калинине взбрыки при полном погружении УПЗ (не отключал централиные приводы)
            else:
                hgrp[0] = 0.0
            SetHgrp3(hgrp, NH, Hnums)
            time_stop = 3600 * 4
            print('upz')
            oldtime = 0
            v.YMFAST = 800
            v.step()
            while float(v.YMTIME) < time_stop:
                if float(v.YMTIME) > 2400 and float(v.YMTIME) > oldtime + float(v.YMFAST) * 0.25:
                    koeff = get_koeffs(bgk, erheat, wnom, plant, Hnums, Icam0)
                    transient.append(koeff)
                    oldtime = float(v.YMTIME)
                v.step()

        hgrp = np.ones(len(Hnums))*1.03
        for h_num in reversed(range(1, len(Hnums))):
            h_max = 6 if h_num == len(Hnums)-1 else 5
            for h_step in reversed(range(h_max)):
                print('SUZ {}'.format(h_step))

                suz_step = (len(Hnums) - h_num - 1) * 5 + (5 - h_step)
                if float(v.YMBOR_COR) > 0.15 * suz_step:
                    set_state_for_dyn(w_cur, NH, Hnums)
                    gs.fix()

                    v.YMFAST = 1
                    v.step()

                    hgrp[h_num] = h_step / 5.
                    SetHgrp3(hgrp, NH, Hnums)
                    v.step()

                    v["#YM#YMTIME"] = 0

                    time_stop = 3600 * 2.
                    if h_num == len(Hnums)-1 and h_step in [0, 3]:
                        time_stop = 3600 * 24.

                    v.YMFAST = 800
                    v.step()

                    print('START')
                    oldtime = 0

                    while float(v.YMTIME) < time_stop:
                        if float(v.YMBOR_COR) > 0.1:
                            if float(v.YMTIME) > 1800 and float(v.YMTIME) > oldtime + float(v.YMFAST) * 0.25:
                                print(burn, w_cur, h_num, h_step, float(v.YMTIME), float(v.YMBOR_COR))
                                koeff = get_koeffs(bgk, erheat, wnom, plant, Hnums, Icam0)
                                transient.append(koeff)
                                oldtime = float(v.YMTIME)
                        else:
                            break
                        v.step()
                    gs.restore()

    LIcam = np.array([point['I'] for point in transient])
    LHgrp = np.array([point['h'] for point in transient])
    LXe = np.array([point['Xe'][: 10] for point in transient]) / xe_norm
    LW = np.array([point['w'][: 10] for point in transient])
    wt = np.array([point['wt'] for point in transient])
    offset = np.array([point['offset'] for point in transient])

    res = {'LIcam': LIcam, 'LHgrp': LHgrp, 'LXe': LXe, 'LW': LW, 'wt': wt, 'offset': offset}
    gs_init.restore()
    v.syncro(0)
    return res


def make_dyn_koef(ncompl, nchan, hdf_in, dims, defaults, burn, use_xe):
    #теперь делаем fit как в математике
    #сначала изготовим функции

    erheat = defaults['erheat'][-1]
    Wnom = defaults['Wnom'][-1]
    nburn = dims['NBurn']

    static_data = {}

    Nw = dims['Nw']
    nfactors = dims['Nv']

    if use_xe:
        polyn = [[0], [1], [2], (1,3), (1,4), (1,5), (1,3,3), (1,4,4), (1,5,5), (1,3,3,3), (1,4,4,4), (1,5,5,5),
                 (1,3,3,3,3), (1,4,4,4,4), (1,5,5,5,5), (0,2), (1,6), (1,7), (1,8), (1,9), (1,10)]
    else:
        polyn = [[0], [1], [2], (1,3), (1,4), (1,5), (1,3,3), (1,4,4), (1,5,5), (1,3,3,3), (1,4,4,4), (1,5,5,5),
                 (1,3,3,3,3), (1,4,4,4,4), (1,5,5,5,5), (0,2)]

    koeFF = ['I1', 'I2', 'I3', 'I2_H4', 'I2_H7', 'I2_H10', 'I2_H4_H4', 'I2_H7_H7', 'I2_H10_H10', 'I2_H4_H4_H4', 'I2_H7_H7_H7', 'I2_H10_H10_H10',
                 'I2_H4_H4_H4_H4', 'I2_H7_H7_H7_H7', 'I2_H10_H10_H10_H10', 'I1_I3', 'I1_Xe1', 'I1_Xe2', 'I1_Xe3', 'I1_Xe4', 'I1_Xe5']

    b_line = np.linspace(0, float(burn), nburn)
    for hdf, b in zip(hdf_in, b_line):
        print(hdf, b)
        result = [[{"v2w": np.zeros((Nw, nfactors))} for ichan in range(nchan)] for icompl in range(ncompl)]
        with h5.File(hdf, 'r') as f5:
            LIcam = f5['LIcam'][:]
            LHgrp = f5['LHgrp'][:]
            LXe = f5['LXe'][:]
            LW = f5['LW'][:]

        for icompl, ichan in itertools.product(range(ncompl), range(nchan)):
        #for iii in range(1):
        #    icompl = 0
        #    ichan = 0
            ftdat = np.hstack((LIcam[:, icompl, ichan], LHgrp, LXe))

            test_polyn = []
            for pol in polyn:
                v_el = np.ones(ftdat.shape[0])
                for ind in pol:
                    v_el *= ftdat[:, ind]
                test_polyn.append(v_el)
                

            # lasso = Lasso(alpha=1e-2, max_iter=5000, fit_intercept=False)
            lasso = Lasso(alpha=1e-2, max_iter=10000, fit_intercept=False, tol=7e-3)
            print('aaa')
            print(test_polyn)
            print('bbb')
            print(np.array(test_polyn).T)
            print('ccc')
            print(LW[:, :Nw])
            lasso.fit(np.array(test_polyn).T, LW[:, :Nw])
            print('ddd')
            print(lasso)
            print('end')
            for iw in range(Nw):
                #print iw, icompl, ichan, b
                fit_koef = list(lasso.coef_[iw])
                print(fit_koef)
                if not use_xe:
                    fit_koef.extend([0]*5)
                result[icompl][ichan]["v2w"][iw, :] = np.array(fit_koef) / Wnom * erheat
            print('here')
            print(result[icompl][ichan]["v2w"][:])
            for i in range(21):
                print(koeFF[i], (result[icompl][ichan]["v2w"][:, i]))

            #for icompl, ichan in itertools.product(range(ncompl), range(nchan)):
            #    result[icompl][ichan] = result[0][0]
        static_data[b] = result
    return static_data

if __name__ == '__main__':
    import glob
    make_dyn_koef(2, 4, glob.glob(r'E:\projects\example\dynamic_5.h5'),
                  ld(r'E:\projects\example\basedims.json'),
                  ld(r'E:\projects\example\defaults.pkl'), 200.0, 1)
    #make_dyn_model(ld(r'F:\Documents\Projects\akpm1200\data\data\novo2\b01\k01\out\dynamic.pkl'),
    #              ld(r'F:\Documents\Projects\akpm1200\data\data\novo2\b01\k01\out\basedims.json'),
    #              ld(r'F:\Documents\Projects\akpm1200\data\data\novo2\b01\k01\out\defaults.pkl'),
    #              200.0, ld(r'F:\Documents\Projects\akpm1200\data\data\novo2\b01\k01\out\norm_koeffs.pkl'),
    #              'novo2', 1, 1)
    

    #parser = argparse.ArgumentParser(prog="generate dynamic data")
    #parser.add_argument("--ncompl", default="2", help=u"number of complects")
    #parser.add_argument("--nchan", default="3", help=u"number of chanals in complect")
    #parser.add_argument("--out", help=u"result")
    #parser.add_argument("--bgk", help=u"bgk", default=[])
    #parser.add_argument("--dims", help=u"dims")
    #parser.add_argument("--defaults", help=u"defaults")
    #parser.add_argument("--initpoint", help=u"initpoint")
    #parser.add_argument("--burn", help=u"burn")
    #parser.add_argument("--plant", help=u"plant")
    #parser.add_argument("--block", help=u"block")
    #parser.add_argument("--kamp", help=u"kamp")
    #parser.add_argument("--hdf_in", help=u"hdf_in", default=[], nargs='*')
    #parser.add_argument("--use_xe", help=u"use_xe", default=1)
    #args = parser.parse_args()
    #ext = os.path.splitext(args.out)[-1]
    #if ext == '.pkl':
    #    res = make_dyn_koef(int(args.ncompl), int(args.nchan), args.hdf_in, ld(args.dims), ld(args.defaults),
    #                        float(args.burn), int(args.use_xe))
    #    sv(res, args.out)
    #elif ext == '.h5':
    #    #res = make_dyn_model(ld(args.bgk), ld(args.dims), ld(args.defaults), ld(args.initpoint), float(args.burn),
    #    #                     ld(args.norm), args.plant, int(args.block), int(args.kamp))
    #    res = make_dyn_model(ld(args.bgk), args.hdf_in)
    #    dump(res, args.out)


