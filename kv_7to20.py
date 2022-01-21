import pickle as pkl
import numpy as np
from mc2py.zgiw import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as us
from mc2py.updn import updn
from mc2py.showkart import *

def coor(n):
	return np.linspace(0,1-1/n,n)+1/n/2

with open(r'Kv/kv.pkl', "rb") as resfile:
    data = pkl.load(resfile)
data = data.reshape(36,7,163)

#kshow(np.sum(data[32], axis=0),show=1)
data = data[:,:,updn - 1]
#kshow(np.sum(data[32], axis=0),show=1)

data_new=[]
for i in range(36):
	data_new.append(data[i].transpose())
data_new = np.array(data_new) 

h = coor(7)
h_ = coor(20)
kv = []
for i in range(36):
	for j in range(163):
		kv_spl = us(h, data_new[i][j], k=3, s=0)
		for q in h_:
			kv.append(float(kv_spl(q)))
kv = np.array(kv).reshape(36,163,20)

with open(r'out_kv/kv_20.pkl', "wb") as resfile:
    pkl.dump(kv, resfile)

for i in range(1):
	for j in range(163):
		plt.plot(h_,kv[i][j],marker='o')
		plt.plot(h,data_new[i][j],marker='*')
		plt.grid(True)
plt.show()