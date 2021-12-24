import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import curve_fit

#E:\projects\Novo2_Zone\MODEL_\MODEL\STATE\st4akpm\novo7_k3_4akpm_static.sta
def approk(x,k,b):
	return k*x+b

def plot_csv_amperage(a):
	a = str(a)
	data = pd.read_csv(a, sep="\t", index_col = False)
	keys = data.keys()
	name = a.split('.')
	name = name[-2]
	name = name.split('\\')
	name = name[-1]
	plt.rcParams.update({'font.size': 14})
	parametr_1, parametr_cov_1 = curve_fit(approk, np.array(data[keys[3]]), np.array(data[keys[1]]))
	parametr_3, parametr_cov_3 = curve_fit(approk, np.array(data[keys[3]]), np.array(data[keys[2]]))
	func_1 = parametr_1[0]*np.array(data[keys[3]])+parametr_1[1]
	func_3 = parametr_3[0]*np.array(data[keys[3]])+parametr_3[1]
	fig = figure(figsize=(40, 20))
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(data[keys[3]]), np.array(data[keys[1]]), 'o', color = np.random.rand(3)) 
	ax1.plot(np.array(data[keys[3]]), func_1, color = np.random.rand(3)) 
	ax1.plot(np.array(data[keys[3]]), np.array(data[keys[2]]), 'o', color = np.random.rand(3)) 
	ax1.plot(np.array(data[keys[3]]), func_3, color = np.random.rand(3))
	ax1.set_xlabel(keys[3] + ', %')
	ax1.set_ylabel(u'$I$')
	ax1.legend((keys[1], 'func_for_' + keys[1], keys[2], 'func_for_' + keys[2]), loc='best')
	ax1.set_title(name + '.csv')
	savefig(r'out\csv_amperage\\' + name + '.png', orientation='landscape', papertype='a4', format='png')
	print()
	print(f'k_{keys[1]} = ', parametr_1[0])
	print(f'k_{keys[2]} = ', parametr_3[0])
	print()
	return data

def plot_csv_amperage_all(a):
	plt.rcParams.update({'font.size': 14})
	fig = figure(figsize=(40, 20))
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	for data in a:
		keys = data.keys()
		parametr_1, parametr_cov_1 = curve_fit(approk, np.array(data[keys[3]]), np.array(data[keys[1]]))
		parametr_3, parametr_cov_3 = curve_fit(approk, np.array(data[keys[3]]), np.array(data[keys[2]]))
		func_1 = parametr_1[0]*np.array(data[keys[3]])+parametr_1[1]
		func_3 = parametr_3[0]*np.array(data[keys[3]])+parametr_3[1]
		ax1.plot(np.array(data[keys[3]]), np.array(data[keys[1]]), 'o', color = np.random.rand(3), label = keys[1])
		ax1.plot(np.array(data[keys[3]]), func_1, color = np.random.rand(3), label = 'func_for_' + keys[1])
		ax1.plot(np.array(data[keys[3]]), np.array(data[keys[2]]), 'o', color = np.random.rand(3), label = keys[2])
		ax1.plot(np.array(data[keys[3]]), func_3, color = np.random.rand(3), label = 'func_for_' + keys[2])
	ax1.legend()
	ax1.set_xlabel(u'offset, %')
	ax1.set_ylabel(u'$I$')
	ax1.set_title(u'all') 
	savefig(r'out\csv_amperage\all.png', orientation='landscape', papertype='a4', format='png')

def plot_csv_half_sum(a):
	a = str(a)
	data = pd.read_csv(a, sep="\t", index_col = False)
	keys = data.keys()
	name = a.split('.')
	name = name[-2]
	name = name.split('\\')
	name = name[-1]
	plt.rcParams.update({'font.size': 14})
	parametr_1, parametr_cov_1 = curve_fit(approk, np.array(data[keys[3]]), (np.array(data[keys[1]]) + np.array(data[keys[2]]))/2)
	func = parametr_1[0]*np.array(data[keys[3]])+parametr_1[1]
	fig = figure(figsize=(40, 20))
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(data[keys[3]]), (np.array(data[keys[1]]) + np.array(data[keys[2]]))/2, 'o', color = np.random.rand(3)) 
	ax1.plot(np.array(data[keys[3]]), func, color = np.random.rand(3)) 
	ax1.set_xlabel(keys[3] + ', %')
	ax1.set_ylabel(u'$(I_{1} + I_{2})/2$')
	ax1.legend((f'({keys[1]} + {keys[2]})/2', 'func'), loc='best')
	ax1.set_title(name + '.csv') 
	savefig(r'out\csv_half_sum\\' + name + '.png', orientation='landscape', papertype='a4', format='png')
	print()
	print(f'k_({keys[1]} + {keys[2]})/2 = ', parametr_1[0])
	print()
	return data

def plot_csv_half_sum_all(a):
	plt.rcParams.update({'font.size': 14})
	fig = figure(figsize=(40, 20))
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	for data in a:
		keys = data.keys()
		parametr_1, parametr_cov_1 = curve_fit(approk, np.array(data[keys[3]]), (np.array(data[keys[1]]) + np.array(data[keys[2]]))/2)
		func = parametr_1[0]*np.array(data[keys[3]])+parametr_1[1]
		ax1.plot(np.array(data[keys[3]]), (np.array(data[keys[1]]) + np.array(data[keys[2]]))/2, 'o', color = np.random.rand(3), label = f'({keys[1]} + {keys[2]})/2') 
		ax1.plot(np.array(data[keys[3]]), func, color = np.random.rand(3), label = 'func_for_' + f'({keys[1]} + {keys[2]})/2') 
	ax1.set_xlabel(u'offset, %')
	ax1.set_ylabel(u'$(I_{1} + I_{2})/2$')
	ax1.legend()	
	ax1.set_title('all') 
	savefig(r'out\csv_half_sum\all.png', orientation='landscape', papertype='a4', format='png')
	return data





a = []
b = []
pathdata = Path(r'data')
for i in pathdata.glob('*.csv'):
	a.append(plot_csv_amperage(i))	
	b.append(plot_csv_half_sum(i))
plot_csv_amperage_all(a)
plot_csv_half_sum_all(b)