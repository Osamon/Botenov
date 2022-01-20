import shutil
import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from example.koeff_polynom_test import get_data_model, w_from_v2w_burn_data, get_set_index 

def process(dict_index):
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
    shutil.copyfile(file_out, file_out_model)
    v['YMBLOCKNUMBER_TO_LOAD'] = 2
    v["YM_N_KAMP_TO_LOAD"] = 3
    v["YM_XIPI_LDBRNBEG"] = 1  # Загружаем начало выбранной кампании
    v.step(n=100)
    v['#YM#YMFLAGSTAT'] = 1
    #v["yztcam_key"] += 0.1
    v['ym_flag_ake'] = 0 # Отключение непрерывной тарировки
    v.step()
    w = []
    i1_list = []
    i2_list = []
    i3_list = []
    waknp = []
    wakpm = []
    ymintpow =[]
    burn_list = []
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
        waknp.append((point["i1"]+point["i3"])/2)
        wakpm.append(point["wakpm"])
        ymintpow.append(point["ymintpow"])
        
        w.append(w_from_v2w_burn_data(point, file_out))

        burn = point["burn"]
        print(burn)
        v.step()
    v["YMFLAG_BRN"] = 0
    v["#YM#YMFLAG_BRN"] = 0
    v['ymfast'] = 1
    v.step()
    v.syncro(0)
    
    w = np.array(w)*100
    ymintpow = np.array(ymintpow)
    wakpm = np.array(wakpm)
    waknp = np.array(waknp)

    w_norm = w*(ymintpow[50]/w[50])
    wakpm_norm = wakpm*(ymintpow[50]/wakpm[50])
    waknp_norm = waknp*(ymintpow[50]/waknp[50])

    plt.plot(burn_list[0:],w_norm[0:],label='w')
    plt.plot(burn_list[0:],ymintpow[0:],label='ymintpow')
    plt.plot(burn_list[0:],wakpm_norm[0:],label='wakpm')
    plt.plot(burn_list[0:],waknp_norm[0:],label='waknp')
    plt.legend()
    #plt.show()
    plt.savefig(f"out/power__I1_{dict_index['I1']}_I2_{dict_index['I2']}_I3_{dict_index['I3']}.png")
    plt.clf()
    plt.plot(burn_list[0:],(np.array(i1_list))[0:],label='i1')
    plt.plot(burn_list[0:],(np.array(i2_list))[0:],label='i2')
    plt.plot(burn_list[0:],(np.array(i3_list))[0:],label='i3')
    plt.legend()
    #plt.show()
    plt.savefig(f"out/amperage__I1_{dict_index['I1']}_I2_{dict_index['I2']}_I3_{dict_index['I3']}.png")
    plt.clf()
    return dict_index['I1']

file_in = r'E:/projects/example/akpm_for_simulator_refs.h5' 
file_out = r'E:/projects/example/akpm_for_simulator.h5'
file_out_model = r'E:/projects/NOVO2_Karu/Model_/MODEL/RO_PART/CORE/CPPLIB/INPUT_AKNP/XIPI1/novo2/B02/K03/akpm_for_simulator.h5'

dict_index_0 = {"I1":np.linspace(0,11,100),
                "I2":np.linspace(0,11,100),
                "I3":np.linspace(0,11,100),
}

dict_index = {}
for i in range(len(dict_index_0['I1'])):
    for k,v in dict_index_0.items():
        dict_index[k] = v[i]
    #print(dict_index)
    get_set_index(dict_index, file_in, file_out)
    print('I1 = ', process(dict_index))