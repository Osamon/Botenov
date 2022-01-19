'''
Этот скрипт восстанавливает энерговыделение в активной зоне по v2w_burn_data
А так же формирует файл koeff.csv, в котором содержатся индексы коэффициентов при каждом из слагаемых в полиноме для каждой гармоники при выгорании burn, с которыми они входят в массив v2w_burn_data     
'''
import numpy as np
import glob
import h5py as h5
from dynmodel import *
from mc2py.zgiw import *
import matplotlib.pyplot as plt
import pandas as pd

def get_data_model(v,i1_list,i2_list,i3_list,waknp,wakpm,burn_list):
    #v = getv()
    #v.syncro(1)
    #v.step()
    burn = float(v['YMTIME_BRN'])
    i1 = float(v['yq_i_aknp_chan(1,1)'])
    i2 = float(v['yq_i_aknp_chan(1,2)'])
    i3 = float(v['yq_i_aknp_chan(1,3)'])
    h7 = float(v['yshgrp(7)'])
    h11 = float(v['yshgrp(11)'])
    h12 = float(v['yshgrp(12)'])
    
    waknp.append(float(v["YQ_W_AKNP(1,1)"]))
    wakpm.append(float(v["YQ_AKPM_WT(1,1)"]))
    burn_list.append(float(v["YMTIME_BRN"]))



    #v.syncro(0)

    i1_0 = 4.13828e+09 #Ток нижней камеры в номинале
    i2_0 = 4.54993e+09 #Ток средней камеры в номинале
    i3_0 = 3.29128e+09 #Ток верхней камеры в номинале

    i1 = i1/i1_0
    i2 = i2/i2_0
    i3 = i3/i3_0
    i1_list.append(i1)
    i2_list.append(i2)
    i3_list.append(i3)

    params = {"burn": burn,
              "i1": i1,
              "i2": i2,
              "i3": i3,
              "h7": h7,
              "h11": h11,
              "h12": h12,
    }
    return params

def w_from_v2w_burn_data(params):

    # get_data_model() # получение данных для алгоритма из модели
    try:
        burn = params["burn"]
    except Exception:
        burn = 1
    try:
        i1 = params["i1"]
    except Exception:
        i1 = 1
    try:
        i2 = params["i2"]
    except Exception:
        i2 = 1
    try:
        i3 = params["i3"]
    except Exception:
        i3 = 1
    try:
        h7 = params["h7"]
    except Exception:
        h7 = 1
    try:
        h11 = params["h11"]
    except Exception:
        h11 = 1
    try:
        h12 = params["h12"]
    except Exception:
        h12 = 1

    Nt = 10 #количество точек по выгоранию
    Nw = 10 #Размерность базиса по мощности
    Nv = 21 #Размерность вектора факторов
    Nf = Nw*Nv
    with h5.File(r'E:\projects\example\akpm_for_simulator.h5', "r") as rf:
        aaa = rf['i_norm'][:]
        v2w_burn_data = rf['lasso/0/0/v2w_burn_data'][:]
        w2Kz = rf['w2Kz'][:]
    x0 = v2w_burn_data[0]
    h = v2w_burn_data[1]
    v2w_burn_data = v2w_burn_data[2:]
    data = np.reshape(v2w_burn_data, (Nt,Nf,2))

    #Восстановление энерговыделения в активной зоне по v2w_burn_data
    #w =[]
    #burn_list = np.arange(0, 300, 10)
    #for burn in burn_list:
    out = []
    index = int((burn - x0)/h)
    if index < 0:
        index = 0
    else:
        if index >= Nt - 1:
            index = Nt - 2
    x = ((burn-x0)-h*index)/h
    for i in range(Nf):
        f0  = data[index][i][0];
        f1  = data[index+1][i][0];
        df0 = data[index][i][1]*h;
        df1 = data[index+1][i][1]*h;
        df = f0 - f1;
        out.append(f0 + x*(df0 - x*(3*df + (2*df0 + df1) - x*(2*df + df0 + df1))))
    out = np.reshape(out, (Nw,Nv))
    #koef = make_dyn_koef(2, 4, glob.glob(r'E:\projects\example\dynamic_5.h5'), ld(r'E:\projects\example\basedims.json'), ld(r'E:\projects\example\defaults.pkl'), 0.0, 1)
    #koef = koef[0][0][0]['v2w']
    q = [i1, i2, i3, i2*h7, i2*h11, i2*h12, i2*h7*h7, i2*h11*h11, i2*h12*h12, i2*h7*h7*h7, i2*h11*h11*h11, i2*h12*h12*h12, i2*h7*h7*h7*h7, i2*h11*h11*h11*h11, i2*h12*h12*h12*h12, i1*i3,0,0,0,0,0]
    koeFF = ['I1', 'I2', 'I3', 'I2_H7', 'I2_H11', 'I2_H12', 'I2_H7_H7', 'I2_H11_H11', 'I2_H12_H12', 'I2_H7_H7_H7', 'I2_H11_H11_H11', 'I2_H12_H12_H12',
                 'I2_H7_H7_H7_H7', 'I2_H11_H11_H11_H11', 'I2_H12_H12_H12_H12', 'I1_I3', 'I1_Xe1', 'I1_Xe2', 'I1_Xe3', 'I1_Xe4', 'I1_Xe5']
    koeffs_garm = np.matmul(out, q)
    #w.append(np.sum(np.matmul(w2Kz, koeffs_garm)))
    w = np.sum(np.matmul(w2Kz, koeffs_garm))
    return w
    #plt.plot(burn_list, w)
    #plt.show()
    '''
    #Формирование файла koeff.csv
    itog = {}
    for i in range(0,42,2):
    	asd = []
    	for j in range(i,4200,42):
    		asd.append(j+2)
    	print(koeFF[int(i/2)], '', asd)
    	itog[f'{koeFF[int(i/2)]}'] = asd

    DF = pd.DataFrame.from_records(itog)
    DF.to_csv(r'koeff.csv', sep="\t")
    '''

#w_from_v2w_burn_data()