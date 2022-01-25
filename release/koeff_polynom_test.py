'''
Этот скрипт восстанавливает энерговыделение в активной зоне по v2w_burn_data
А так же формирует файл koeff.csv, в котором содержатся индексы коэффициентов при каждом из слагаемых в полиноме для каждой гармоники при выгорании burn, с которыми они входят в массив v2w_burn_data.
Все функции, непосредственно взаимодействующие с моделью, были протестированы на этой модели: E:/projects/NOVO2_Karu/Model_/MODEL/karru_nv_2giw_execute.bat
'''

import pickle as pkl
import numpy as np
import glob
import h5py as h5
from mc2py.zgiw import *
import matplotlib.pyplot as plt
import pandas as pd
import shutil

'''
koeFF - этот массив содержит названия коэффициентов, стоящих по порядку перед слагаемыми в полиноме мощности, названия этих коэффициентов соответсвуют тем входным сигналам, перед которыми они стоят в полиноме
'''
koeFF = ['I1', 'I2', 'I3', 'I2_H7', 'I2_H11', 'I2_H12', 'I2_H7_H7', 'I2_H11_H11', 'I2_H12_H12', 'I2_H7_H7_H7', 'I2_H11_H11_H11', 'I2_H12_H12_H12',
         'I2_H7_H7_H7_H7', 'I2_H11_H11_H11_H11', 'I2_H12_H12_H12_H12', 'I1_I3', 'I1_Xe1', 'I1_Xe2', 'I1_Xe3', 'I1_Xe4', 'I1_Xe5']

def get_data_model(v):
    '''
    Эта функция позволяет набирать входные данные для работы АКПМ (w_from_v2w_burn_data) из модели.
    Всё работает хорошо, но важно проверять работоспособность используемых ККС для конкретной модели и не забывать изменять значения токов в номинале (i1_0, i2_0, i3_0) для соответствующей модели, блока и кампании
    '''
    #v = getv()
    #v.syncro(1)
    #v.step()
    burn = float(v['YMTIME_BRN'])
    waknp_model = float(v["YQ_W_AKNP(1,1)"])
    wakpm_model = float(v["YQ_AKPM_WT(1,1)"])
    ymintpow = float(v["ymintpow"])
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

    waknp = (i1+i3)/2

    point = {"burn": burn,
              "waknp": waknp,
              "ymintpow": ymintpow,
              "waknp_model": waknp_model,
              "wakpm_model": wakpm_model,
              "i1": i1,
              "i2": i2,
              "i3": i3,
              "h7": h7,
              "h11": h11,
              "h12": h12,
    }
    return point

def w_from_v2w_burn_data(params=None, file_out=r'E:/projects/example/akpm_for_simulator.h5'):
    '''Восстановление энерговыделения в активной зоне по v2w_burn_data (имитация АКПМ),
    всё в этой функции работает хорошо'''

    # get_data_model() # получение данных для алгоритма из модели
    if params != None:
        try:
            burn = params["burn_list"]
        except Exception:
            burn = 1
        try:
            i1 = params["i1_list"]
        except Exception:
            i1 = 1
        try:
            i2 = params["i2_list"]
        except Exception:
            i2 = 1
        try:
            i3 = params["i3_list"]
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
    else:
        burn = 1
        i1 = 1
        i2 = 1
        i3 = 1
        h7 = 1
        h11 = 1
        h12 = 1

    Nt = 10 #количество точек по выгоранию
    Nw = 10 #Размерность базиса по мощности
    Nv = 21 #Размерность вектора факторов
    Nf = Nw*Nv

    with h5.File(file_out, "r") as rf:
        v2w_burn_data = rf['lasso/0/0/v2w_burn_data'][:]
        w2Kz = rf['w2Kz'][:]
    x0 = v2w_burn_data[0]
    h = v2w_burn_data[1]
    v2w_burn_data = v2w_burn_data[2:]
    data = np.reshape(v2w_burn_data, (Nt,Nf,2))

    #wakpm =[]
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
    #wakpm.append(np.sum(np.matmul(w2Kz, koeffs_garm)))
    #plt.plot(burn_list, wakpm)
    #plt.show()
    wakpm = np.sum(np.matmul(w2Kz, koeffs_garm))
    return wakpm

def get_set_index(dict_index=None, file_in=r'E:/projects/example/akpm_for_simulator_refs.h5', file_out=r'E:/projects/example/akpm_for_simulator.h5'):
    '''Формирование файла koeff.csv'''
    '''Эта функция позволяет менять коэффициенты при любых слагаемых в полиноме мощности в каждой из 10 гармоник, во всех десяти точках по выгоранию.
    На выходе получаем изменённый файл akpm_for_simulator_refs.h5. Всё рабоатет хорошо'''

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
                v2w_burn_data[itog[k][i]+1] = v2w_burn_data[itog[k][i]+1]*v

        with h5.File(file_out, 'a') as f:
            dset = f['lasso/0/0/v2w_burn_data']
            dset[:] = v2w_burn_data

def process():
    '''Эта функция просто набирает данные из модели, в процессе выгорание активной зоны.
    Всё работает хорошо, но для корректной работы необходимо проверять актуальнсоть ККС'''

    print('Работаю с моделью')
    v=getv()
    v.load_state(r"KIRIN/B01_K01_Nominal_BOC") #B02_K02_Nominal_BOC_0-01eff_days.sta")
    v.step(n=1)
    v.syncro(1)
    v["KEY1_SHUNT_AZ1"] += 1
    v["KEY_SHUNT_URB"] += 1
    v["KEY_SHUNT_PZ1"] += 1
    v["KEY_SHUNT_PZ2"] += 1
    v["KLARM_POS"] = 1
    if float(v['YZTABLINKNKS']) == False:
        v["YZKEYLINK2"] += 0.2 # переходим в режим НКС
    v['#YM#YMFLAGSTAT'] = 1
    v.step(n=5)
    v.YZBORMODE = 2
    v.step()
    #shutil.copyfile(file_out, file_out_model)
    v['YMBLOCKNUMBER_TO_LOAD'] = 2
    v["YM_N_KAMP_TO_LOAD"] = 2
    v["YM_XIPI_LDBRNBEG"] = 1  # Загружаем начало выбранной кампании
    v.step(n=100)
    v['#YM#YMFLAGSTAT'] = 1
    #v["yztcam_key"] += 0.1
    v['ym_flag_ake'] = 0 # Отключение непрерывной тарировки
    v.step()
    burn_list = []
    i1_list = []
    i2_list = []
    i3_list = []
    waknp = []
    waknp_model = []
    wakpm_model = []
    ymintpow =[]
    h7 = []
    h11 = []
    h12 = []

    v['ymfast'] = 50000
    burn = float(v['YMTIME_BRN'])
    v["YMFLAG_BRN"] = 1
    v["#YM#YMFLAG_BRN"] = 1
    v.step()
    while burn < 100:
        point = get_data_model(v)

        burn_list.append(point["burn"])
        i1_list.append(point["i1"])
        i2_list.append(point["i2"])
        i3_list.append(point["i2"])
        waknp.append(point["waknp"])
        waknp_model.append(point["waknp_model"])
        wakpm_model.append(point["wakpm_model"])
        ymintpow.append(point["ymintpow"])
        h7.append(point["h7"])
        h11.append(point["h11"])
        h12.append(point["h12"])

        burn = point["burn"]
        print(burn)
        v.step()

    v["YMFLAG_BRN"] = 0
    v["#YM#YMFLAG_BRN"] = 0
    v['ymfast'] = 1
    v.step()
    v.syncro(0)

    dict_model = {"burn_list":burn_list,
                  "i1_list":i1_list,
                  "i2_list":i2_list,
                  "i3_list":i3_list,
                  "waknp":waknp,
                  "waknp_model":waknp_model,
                  "wakpm_model":wakpm_model,
                  "ymintpow":ymintpow,
                  "h7": h7,
                  "h11": h11,
                  "h12": h12,
    }

    with open(r'out/pkl/from_model.pkl', "wb") as resfile:
        pkl.dump(dict_model, resfile)
    print('Модель отработала')

def process_upz():
    '''Эта функция моделирует срабатывание УПЗ.
    Всё работает хорошо, но для корректной работы необходимо проверять актуальнсоть ККС'''

    print('Работаю с моделью')
    v=getv()
    v.load_state(r"KIRIN/B01_K01_Nominal_BOC") #B02_K02_Nominal_BOC_0-01eff_days.sta")
    v.syncro(1)
    v.step(n=1)

    v.step()
    burn_list = []
    i1_list = []
    i2_list = []
    i3_list = []
    waknp = []
    waknp_model = []
    wakpm_model = []
    ymintpow = []
    h7 = []
    h11 = []
    h12 = []
    step = 0
    v.step()
    while step < 10:
        point = get_data_model(v)

        burn_list.append(point["burn"])
        i1_list.append(point["i1"])
        i2_list.append(point["i2"])
        i3_list.append(point["i2"])
        waknp.append(point["waknp"])
        waknp_model.append(point["waknp_model"])
        wakpm_model.append(point["wakpm_model"])
        ymintpow.append(point["ymintpow"])
        h7.append(point["h7"])
        h11.append(point["h11"])
        h12.append(point["h12"])

        print(step)
        step+=1
        v.step()

    v['11JDY01CH403_SBV1'] = 1
    print('Сработала УПЗ')
    step = 0
    v.step()
    while step < 100:
        point = get_data_model(v)

        burn_list.append(point["burn"])
        i1_list.append(point["i1"])
        i2_list.append(point["i2"])
        i3_list.append(point["i2"])
        waknp.append(point["waknp"])
        waknp_model.append(point["waknp_model"])
        wakpm_model.append(point["wakpm_model"])
        ymintpow.append(point["ymintpow"])
        h7.append(point["h7"])
        h11.append(point["h11"])
        h12.append(point["h12"])

        print(step)
        step+=1
        v.step()

    v.step()
    v.syncro(0)

    dict_model = {"burn_list":burn_list,
                  "i1_list":i1_list,
                  "i2_list":i2_list,
                  "i3_list":i3_list,
                  "waknp":waknp,
                  "waknp_model":waknp_model,
                  "wakpm_model":wakpm_model,
                  "ymintpow":ymintpow,
                  "h7": h7,
                  "h11": h11,
                  "h12": h12,
    }

    with open(r'out/upz/pkl/from_model.pkl', "wb") as resfile:
        pkl.dump(dict_model, resfile)
    print('Модель отработала')

def plot(dict_index, wakpm):
    '''Эта функция нужна просто для отрисовких данных.
    Всё работает хорошо'''

    with open(r'out/upz/pkl/from_model.pkl', "rb") as resfile:
        dict_model = pkl.load(resfile)
    burn_list = dict_model["burn_list"]
    i1_list = dict_model["i1_list"]
    i2_list = dict_model["i2_list"]
    i3_list = dict_model["i3_list"]
    waknp = dict_model["waknp"]
    waknp_model = dict_model["waknp_model"]
    wakpm_model = dict_model["wakpm_model"]
    ymintpow = dict_model["ymintpow"]

    wakpm = np.array(wakpm)*100
    print(wakpm)

    ymintpow = np.array(ymintpow)
    waknp = np.array(waknp)
    wakpm_model = np.array(wakpm_model)

    wakpm_norm = wakpm*(ymintpow[50]/wakpm[50])
    waknp_norm = waknp*(ymintpow[50]/waknp[50])
    wakpm_model_norm = wakpm_model*(ymintpow[50]/wakpm_model[50])

    step_list = np.linspace(0,110,110)
    begin = 0
    plt.plot(step_list[begin:],wakpm_norm[begin:],label='wakpm')#,marker='o')
    plt.plot(step_list[begin:],ymintpow[begin:],label='ymintpow')
    plt.plot(step_list[begin:],waknp_norm[begin:],label='waknp')
    plt.plot(step_list[begin:],wakpm_model_norm[begin:],label='wakpm_model')
    plt.legend()
    plt.savefig(f"out/upz/power__I2_H7_{dict_index['I2_H7']}.png")
    plt.clf()
    """
    plt.plot(burn_list[begin:],(np.array(i1_list))[begin:],label='i1')
    plt.plot(burn_list[begin:],(np.array(i2_list))[begin:],label='i2')
    plt.plot(burn_list[begin:],(np.array(i3_list))[begin:],label='i3')
    plt.legend()
    plt.savefig(f"out/amperage__I1_{dict_index['I1']}_I2_{dict_index['I2']}_I3_{dict_index['I3']}.png")
    plt.clf()
    """
def sko(ymintpow, wakpm):
    '''Эта функция считает СКО, всё работает хорошо'''

    return np.sqrt(np.sum((np.array(ymintpow) - np.array(wakpm)*100)**2 ) / (len(ymintpow) - 1))

def dependence(dict_index_0, key, skolist):
    '''Это для выявления зависимости СКО от разных h5.
    Здесь просто отрисовывается множество графиков зависимости СКО от изменения коэффициентов при различных слагаемых в полиноме мощности.
    Всё работает хорошо'''
    fig = plt.figure()
    ax = fig.add_subplot(111)#, projection='3d')
    ax.plot(dict_index_0[key], skolist, "bo")
    ax.legend('sko')
    ax.set_xlabel(f'Коэф при {key}')
    #ax.set_ylabel(u'Коэф при I1')
    ax.set_ylabel(u'СКО')
    #plt.plot(dict_index_0["I1"], skolist,label='sko', marker='o')
    plt.savefig(f"out/upz/pkl/SKO_{key}.png")
    print('Картинка СКО готова')
    #plt.show()
    #plt.clf()

#print(__doc__)
#print(get_set_index.__doc__)
#print(w_from_v2w_burn_data.__doc__)
#print(plot.__doc__)

#dict_index = {}
#for i in koeFF:
#    dict_index[i] = 1
#dict_index = {'I3':1}
#get_set_index(dict_index)
#print(dict_index)
#print(w_from_v2w_burn_data())
#for i in koeFF:
#    dict_index[i] = 2
#dict_index = {'I3':10}
#get_set_index(dict_index)
#print(dict_index)
#print(w_from_v2w_burn_data())