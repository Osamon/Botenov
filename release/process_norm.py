import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from example.koeff_polynom_test import *

file_in = r'E:/projects/example/akpm_for_simulator_refs.h5' 
file_out = r'E:/projects/example/akpm_for_simulator.h5'
file_out_model = r'E:/projects/NOVO2_Karu/Model_/MODEL/RO_PART/CORE/CPPLIB/INPUT_AKNP/XIPI1/novo2/B02/K03/akpm_for_simulator.h5'

dict_index_0 = {"I1":np.linspace(0,9.9,10),
                "I2":np.linspace(0,9.9,10),
                "I3":np.linspace(0,9.9,10),
}

dict_index = {}
params = {}
skolist = []
wakpm = []
#process()
print('Сохраняю картинки')
with open(r'out/pkl/from_model.pkl', "rb") as resfile:
    dict_model = pkl.load(resfile)

for i in range(len(dict_index_0['I1'])):
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

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(dict_index_0["I1"], dict_index_0["I3"], skolist, label='sko')

#plt.plot(dict_index_0["I1"], skolist,label='sko', marker='o')
plt.savefig("out/pkl/SKO.png")
print('Картинка СКО готова')
plt.show()
#plt.clf()