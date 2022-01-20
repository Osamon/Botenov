import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from example.koeff_polynom_test import get_data_model, w_from_v2w_burn_data 

v=getv()
v.syncro(1)
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
    
    w.append(w_from_v2w_burn_data(point))

    burn = float(v['YMTIME_BRN'])
    print(burn)
    v.step()
v['ymfast'] = 1
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
plt.savefig('F1.png')
plt.clf()
plt.plot(burn_list[0:],(np.array(i1_list))[0:],label='i1')
plt.plot(burn_list[0:],(np.array(i2_list))[0:],label='i2')
plt.plot(burn_list[0:],(np.array(i3_list))[0:],label='i3')
plt.legend()
#plt.show()
plt.savefig('F2.png')