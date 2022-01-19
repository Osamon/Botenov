import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from koeff_polynom_test import get_data_model, w_from_v2w_burn_data 

v=getv()
v.syncro(1)
v.step()
w = []
i1_list = []
i2_list = []
i3_list = []
waknp = []
wakpm = []
burn_list = []
v['ymfast'] = 50000
burn = float(v['YMTIME_BRN'])
v.step()
while burn < 100:
    w.append(w_from_v2w_burn_data(get_data_model(v,i1_list,i2_list,i3_list,waknp,wakpm,burn_list)))
    burn = float(v['YMTIME_BRN'])
    print(burn)
    v.step()
v['ymfast'] = 1
v.syncro(0)
plt.plot(burn_list[1:],(np.array(w)/0.93*100)[1:],label='w')
plt.plot(burn_list[1:],(np.array(wakpm))[1:],label='wakpm')
plt.plot(burn_list[1:],(np.array(waknp))[1:],label='waknp')
plt.legend()
#plt.show()
plt.savefig('F1.png')
plt.clf()
plt.plot(burn_list[1:],(np.array(i1_list))[1:],label='i1')
plt.plot(burn_list[1:],(np.array(i2_list))[1:],label='i2')
plt.plot(burn_list[1:],(np.array(i3_list))[1:],label='i3')
plt.legend()
#plt.show()
plt.savefig('F2.png')