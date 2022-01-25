import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from example.koeff_polynom_test import *

file_in = r'E:/projects/example/akpm_for_simulator_refs.h5' 
file_out = r'E:/projects/example/akpm_for_simulator.h5'
file_out_model = r'E:/projects/NOVO2_Karu/Model_/MODEL/RO_PART/CORE/CPPLIB/INPUT_AKNP/XIPI1/novo2/B02/K03/akpm_for_simulator.h5'

'''
"""Это для выявления зависимости СКО от разных h5"""
value = np.linspace(0,9.9,10)

koeFF = ['I1', 'I2', 'I3', 'I2_H7', 'I2_H11', 'I2_H12', 'I2_H7_H7', 'I2_H11_H11', 'I2_H12_H12', 'I2_H7_H7_H7', 'I2_H11_H11_H11', 'I2_H12_H12_H12',
         'I2_H7_H7_H7_H7', 'I2_H11_H11_H11_H11', 'I2_H12_H12_H12_H12', 'I1_I3', 'I1_Xe1', 'I1_Xe2', 'I1_Xe3', 'I1_Xe4', 'I1_Xe5']

dict_index_0 = {}
#process()
#process_upz()
for p in range(len(koeFF)):
    dict_index_0[f'{koeFF[p]}']=value
    dict_index = {}
    params = {}
    skolist = []
    wakpm = []
    print('Сохраняю картинки')
    with open(r'out/upz/pkl/from_model.pkl', "rb") as resfile:
        dict_model = pkl.load(resfile)

    for i in range(len(dict_index_0[f'{koeFF[p]}'])):
        for k,v in dict_index_0.items():
            dict_index[k] = v[i]
        get_set_index(dict_index, file_in, file_out)
        for j in range(len(dict_model["burn_list"])):
            for kk,vv in dict_model.items():
                params[kk] = vv[j]
            wakpm.append(w_from_v2w_burn_data(params, file_out))
        skolist.append(sko(dict_model["ymintpow"], wakpm))
        #plot(dict_index, wakpm)
        wakpm = []
    dependence(dict_index_0, koeFF[p], skolist)
    dict_index_0 = {}
'''






"""Это для постороения 3d графиков"""
x = np.linspace(0,9.9,10)
y = np.linspace(0,9.9,10)
xv, yv = np.meshgrid(x, y)
xv = xv.reshape(100)
yv = yv.reshape(100)

dict_index_0 = {#"I1":xv,
                #"I2":np.linspace(0,9.9,10),
                "I2_H7":np.linspace(0,9.9,10),
}

dict_index = {}
params = {}
skolist = []
wakpm = []
#process()
#process_upz()
print('Сохраняю картинки')
with open(r'out/upz/pkl/from_model.pkl', "rb") as resfile:
    dict_model = pkl.load(resfile)

for i in range(len(dict_index_0["I2_H7"])):
    for k,v in dict_index_0.items():
        dict_index[k] = v[i]
    get_set_index(dict_index, file_in, file_out)
    for j in range(len(dict_model["burn_list"])):
        for kk,vv in dict_model.items():
            params[kk] = vv[j]
        wakpm.append(w_from_v2w_burn_data(params, file_out))
    skolist.append(sko(dict_model["ymintpow"], wakpm))
    plot(dict_index, wakpm)
    wakpm = []
"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(dict_index_0["I1"], dict_index_0["I3"], skolist, "bo")
ax.legend('sko')
ax.set_xlabel(u'Коэф при I1')
ax.set_ylabel(u'Коэф при I3')
ax.set_zlabel(u'СКО')
#plt.plot(dict_index_0["I1"], skolist,label='sko', marker='o')
plt.savefig(f"out/upz/3d/SKO_I1_I3.png")
print('Картинка СКО готова')
plt.show()
"""