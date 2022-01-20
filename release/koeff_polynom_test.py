'''
Этот скрипт восстанавливает энерговыделение в активной зоне по v2w_burn_data
А так же формирует файл koeff.csv, в котором содержатся индексы коэффициентов при каждом из слагаемых в полиноме для каждой гармоники при выгорании burn, с которыми они входят в массив v2w_burn_data
'''
import numpy as np
import glob
import h5py as h5
from mc2py.zgiw import *
import matplotlib.pyplot as plt
import pandas as pd
import shutil

def get_data_model(v):
    #v = getv()
    #v.syncro(1)
    #v.step()
    burn = float(v['YMTIME_BRN'])
    wakpm = float(v["YQ_AKPM_WT(1,1)"])
    ymintpow = float(v["ymintpow"])
    waknp_model = float(v["YQ_W_AKNP(1,1)"])
    i1 = float(v['yq_i_aknp_chan(1,1)'])
    i2 = float(v['yq_i_aknp_chan(1,2)'])
    i3 = float(v['yq_i_aknp_chan(1,3)'])
    h7 = float(v['yshgrp(7)'])
    h11 = float(v['yshgrp(11)'])
    h12 = float(v['yshgrp(12)'])

    #v.syncro(0)

    i1_0 = 2.65875e+09#4.13828e+09 #Ток нижней камеры в номинале
    i2_0 = 3.02035e+09#4.54993e+09 #Ток средней камеры в номинале
    i3_0 = 2.33893e+09#3.29128e+09 #Ток верхней камеры в номинале

    i1 = i1/i1_0
    i2 = i2/i2_0
    i3 = i3/i3_0

    params = {"burn": burn,
              "wakpm": wakpm,
              "ymintpow": ymintpow,
              "waknp_model": waknp_model,
              "i1": i1,
              "i2": i2,
              "i3": i3,
              "h7": h7,
              "h11": h11,
              "h12": h12,
    }
    return params

def w_from_v2w_burn_data(params, file_out=r'E:/projects/example/akpm_for_simulator.h5'):

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

    with h5.File(file_out, "r") as rf:
        #aaa = rf['i_norm'][:]
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
    koeffs_garm = np.matmul(out, q)
    #w.append(np.sum(np.matmul(w2Kz, koeffs_garm)))
    w = np.sum(np.matmul(w2Kz, koeffs_garm))
    return w
    #plt.plot(burn_list, w)
    #plt.show()
    
def get_set_index(dict_index=None, file_in=r'E:/projects/example/akpm_for_simulator_refs.h5', file_out=r'E:/projects/example/akpm_for_simulator.h5'):
    #Формирование файла koeff.csv
    koeFF = ['I1', 'I2', 'I3', 'I2_H7', 'I2_H11', 'I2_H12', 'I2_H7_H7', 'I2_H11_H11', 'I2_H12_H12', 'I2_H7_H7_H7', 'I2_H11_H11_H11', 'I2_H12_H12_H12',
             'I2_H7_H7_H7_H7', 'I2_H11_H11_H11_H11', 'I2_H12_H12_H12_H12', 'I1_I3', 'I1_Xe1', 'I1_Xe2', 'I1_Xe3', 'I1_Xe4', 'I1_Xe5']
    itog = {}
    for i in range(0,42,2):
    	asd = []
    	for j in range(i,4200,42):
    		asd.append(j+2)
    	if dict_index is None:
            print(koeFF[int(i/2)], '', asd)
    	itog[f'{koeFF[int(i/2)]}'] = asd
    
    DF = pd.DataFrame.from_records(itog)
    DF.to_csv(r'koeff.csv', sep="\t")

    shutil.copyfile(file_in, file_out)
    
    if dict_index != None:
        
        with h5.File(file_in, "r") as rf:
            v2w_burn_data = rf['lasso/0/0/v2w_burn_data'][:]

        for k,v in dict_index.items():
            for i in range(len(itog[k])):
                v2w_burn_data[itog[k][i]] = v2w_burn_data[itog[k][i]]*v

        with h5.File(file_out, 'a') as f:
            dset = f['lasso/0/0/v2w_burn_data']
            dset[:] = v2w_burn_data