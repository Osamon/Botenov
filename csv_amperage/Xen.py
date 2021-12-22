from mc2py.zgiw import *
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import curve_fit

#E:\projects\Novo2_Zone\MODEL_\MODEL\STATE\st4akpm\novo7_k3_4akpm_static.sta
def approk(x,k,b):
	return k*x+b

plot = 1
if plot == False:
	v = getv()
	v.syncro(1)
	v.load_state(r"st4akpm\novo7_k3_4akpm_static.sta")
	v.step(n=20)
	v['#YM#YMFLAGSTAT'] = 1
	v.step(n=1)
	
	end_camp = 1
	
	if end_camp == True:
		v['YM_XIPI_LDBRNEND'] = 1
		v.step(n=1)
		v['ym_flag_ake'] = 0
		v.step(n=1)
		v["YZREGRPINS_KEY"] = True
		v.step(n=1)
		v["#YM#YZREGRPINS_KEY"] = True
		v.step(n=1)
		v["#SETN#YZREGRPINS_KEY"] = True
		v.step(n=1)
		v["#YV#YZREGRPINS_KEY"] = True
		v.step(n=1)


	# Ксеноновые колебания при постоянной мощности
	# for h12end_step in np.linspace(0.84, 1.03, 2):
	#set_state_for_bgk(w_cur, NH, Hnums)
	
	NH = 12
	Hnums = [6, 10, 11]

	

	v.YMFAST = 1
	v.step(n=1)           
	# time_stop = 5*24*3600 * (w_cur/100.) # Ксеноновые колебания ждем несколько суток на номинале. На малой мощности нет смысла их ловить долго
	if end_camp == True:
		time_stop = 7*24*3600
	else:
		time_stop = 5*24*3600
	v.ymflag_xe = 1
	hgrp = np.ones(len(Hnums))*1.03
	hgrp[-1] = 0.75
	SetHgrp3(hgrp, NH, Hnums)
	#print('Xe-procces with H12end=' + str(h12end_step))
	v.step(n=20)
	hgrp = np.ones(len(Hnums))*1.03
	hgrp[-1] = 1.03
	SetHgrp3(hgrp, NH, Hnums)
	oldtime = 0
	v["#YM#YMTIME"] = 0
	v.step(n=1)
	v.YMFAST = 2000

	data_res = {}
	data_res['time'] = []
	data_res['ymintpow'] = []
	data_res['YQ_i_AKNP_chan_11'] = []
	data_res['YQ_i_AKNP_chan_12'] = []
	data_res['YQ_i_AKNP_chan_13'] = []
	data_res['YQ_W_AKNP_11'] = [] # Некорректированная мощность по АКНП
	data_res['YQ_AKPM_WT_11'] = [] # Корректированная мощность по АКНП
	data_res['offset_svrk'] = []
	data_res['offset_akpm'] = []

	while float(v.YMTIME) < time_stop:
		data_res['time'].append(float(v['YMTIME']))
		data_res['ymintpow'].append(float(v['ymintpow']))
		data_res['YQ_i_AKNP_chan_11'].append(float(v['yq_i_aknp_chan(1,1)']))
		data_res['YQ_i_AKNP_chan_12'].append(float(v['yq_i_aknp_chan(1,2)']))
		data_res['YQ_i_AKNP_chan_13'].append(float(v['yq_i_aknp_chan(1,3)']))
		data_res['YQ_W_AKNP_11'].append(float(v['yq_w_aknp(1, 1)']))
		data_res['YQ_AKPM_WT_11'].append(float(v['yq_akpm_wt(1, 1)']))
		data_res['offset_svrk'].append(float(v['YMOFFSETn_']))
		data_res['offset_akpm'].append(float(v['yq_akpm_offset(1, 1)'])*100)
		# bgk.append(wnom,plant=plant,Icam0=Icam0,Hnums=Hnums)
		v.step(n=1)
	v.ymflag_xe = 0
	v.syncro(0)

	begin = 100 
	data_res['time'] = data_res['time'][begin:]
	data_res['ymintpow'] = data_res['ymintpow'][begin:] 
	data_res['YQ_i_AKNP_chan_11'] = data_res['YQ_i_AKNP_chan_11'][begin:]
	data_res['YQ_i_AKNP_chan_12'] = data_res['YQ_i_AKNP_chan_12'][begin:]
	data_res['YQ_i_AKNP_chan_13'] = data_res['YQ_i_AKNP_chan_13'][begin:]
	data_res['YQ_W_AKNP_11'] = data_res['YQ_W_AKNP_11'][begin:]
	data_res['YQ_AKPM_WT_11'] = data_res['YQ_AKPM_WT_11'][begin:]
	data_res['offset_svrk'] = data_res['offset_svrk'][begin:]
	data_res['offset_akpm'] = data_res['offset_akpm'][begin:]

	df = pd.DataFrame(data_res)
	if end_camp == True:
		df.to_csv((r"out\end_camp.csv"), sep="\t", index=False)
		with open(r'out\data_res_end.pkl', "wb") as resfile:
			pickle.dump(data_res, resfile)
		print('Можно рисовать')
	else:
		df.to_csv((r"out\begin_camp.csv"), sep="\t", index=False)
		with open(r'out\data_res_begin.pkl', "wb") as resfile:
			pickle.dump(data_res, resfile)
		print('Можно рисовать') 


csv = True
if csv:
	
	dataIoff_block = pd.read_csv(r'data\dataIoff_block.csv', sep="\t", index_col = False)
	plt.rcParams.update({'font.size': 14})
	parametr_1, parametr_cov_1 = curve_fit(approk, np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I1block']))
	parametr_3, parametr_cov_3 = curve_fit(approk, np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I3block']))
	func_11 = parametr_1[0]*np.array(dataIoff_block['offset_block'])+parametr_1[1]
	func_31 = parametr_3[0]*np.array(dataIoff_block['offset_block'])+parametr_3[1]
	fig = figure(figsize=(40, 20))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I1block']), 'ro', np.array(dataIoff_block['offset_block']), func_11, 'k', 
			 np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I3block']), 'bo', np.array(dataIoff_block['offset_block']), func_31, 'g')
	ax1.set_xlabel(u'offset_block, %')
	ax1.set_ylabel(u'$I$')
	ax1.legend((u'I1block', u'func_11', u'I3block', u'func_31'), loc='best')
	ax1.set_title(u'dataIoff_block.csv') 
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	savefig(r'out\csv\dataIoff_block.png', orientation='landscape', papertype='a4', format='png')
	print()
	print(parametr_1[0])
	print(parametr_3[0])
	print()
	

	dataIoff_model1 = pd.read_csv(r'data\dataIoff_model1.csv', sep="\t", index_col = False)
	plt.rcParams.update({'font.size': 14})
	parametr_1, parametr_1_cov = curve_fit(approk, np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I1model_1']))
	parametr_3, parametr_3_cov = curve_fit(approk, np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I3model_1']))
	func_12 = parametr_1[0]*np.array(dataIoff_model1['offset_model1'])+parametr_1[1]
	func_32 = parametr_3[0]*np.array(dataIoff_model1['offset_model1'])+parametr_3[1]
	fig = figure(figsize=(40, 20))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I1model_1']), 'ro', np.array(dataIoff_model1['offset_model1']), func_12, 'k',
			 np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I3model_1']), 'bo', np.array(dataIoff_model1['offset_model1']), func_32, 'g')
	ax1.set_xlabel(u'offset_block, %')
	ax1.set_ylabel(u'$I$')
	ax1.legend((u'I1model_1', u'func_12', u'I3model_1', u'func_32'), loc='best')
	ax1.set_title(u'dataIoff_model1.csv') 
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	savefig(r'out\csv\dataIoff_model1.png', orientation='landscape', papertype='a4', format='png')
	print()
	print(parametr_1[0])
	print(parametr_3[0])
	







	data_begin_camp = pd.read_csv(r'data\begin_camp.csv', sep="\t", index_col = False)
	plt.rcParams.update({'font.size': 14})
	YQ_i_AKNP_chan_11_cp = np.sum(data_begin_camp['YQ_i_AKNP_chan_11'])/len(data_begin_camp['YQ_i_AKNP_chan_11'])
	YQ_i_AKNP_chan_13_cp = np.sum(data_begin_camp['YQ_i_AKNP_chan_13'])/len(data_begin_camp['YQ_i_AKNP_chan_13'])
	YQ_i_AKNP_chan_11 = np.array(data_begin_camp['YQ_i_AKNP_chan_11'])/YQ_i_AKNP_chan_11_cp
	YQ_i_AKNP_chan_13 = np.array(data_begin_camp['YQ_i_AKNP_chan_13'])/YQ_i_AKNP_chan_13_cp
	parametr_1, parametr_cov_1 = curve_fit(approk, np.array(data_begin_camp['offset_svrk']), YQ_i_AKNP_chan_11)
	parametr_3, parametr_cov_3 = curve_fit(approk, np.array(data_begin_camp['offset_svrk']), YQ_i_AKNP_chan_13)
	func_13 = parametr_1[0]*np.array(data_begin_camp['offset_svrk'])+parametr_1[1]
	func_33 = parametr_3[0]*np.array(data_begin_camp['offset_svrk'])+parametr_3[1]
	fig = figure(figsize=(40, 20))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I1block']), 'ro', np.array(dataIoff_block['offset_block']), func_11, 'k', 
			 np.array(dataIoff_block['offset_block']), np.array(dataIoff_block['I3block']), 'bo', np.array(dataIoff_block['offset_block']), func_31, 'g',
			 np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I1model_1']), 'co', np.array(dataIoff_model1['offset_model1']), func_12, 'k',
			 np.array(dataIoff_model1['offset_model1']), np.array(dataIoff_model1['I3model_1']), 'mo', np.array(dataIoff_model1['offset_model1']), func_32, 'g',
			 np.array(data_begin_camp['offset_svrk']) - 4.5, YQ_i_AKNP_chan_11, 'yo', np.array(data_begin_camp['offset_svrk']) - 4.5, func_13, 'k', 
			 np.array(data_begin_camp['offset_svrk']) - 4.5, YQ_i_AKNP_chan_13, 'go', np.array(data_begin_camp['offset_svrk']) - 4.5, func_33, 'b')
	ax1.set_xlabel(u'offset, %')
	ax1.set_ylabel(u'$I$')
	ax1.legend((u'I1block', u'func_11', u'I3block', u'func_31', u'I1model_1', u'func_12', u'I3model_1', u'func_32', u'YQ_i_AKNP_chan_11', u'func_13', u'YQ_i_AKNP_chan_13', u'func_33'), loc='best')
	ax1.set_title(u'all') 
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	savefig(r'out\csv\all.png', orientation='landscape', papertype='a4', format='png')
	print()
	print(parametr_1[0])
	print(parametr_3[0])






	fig = figure(figsize=(40, 20))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(np.array(data_begin_camp['offset_svrk']), YQ_i_AKNP_chan_11, 'ro', np.array(data_begin_camp['offset_svrk']), func_13, 'k', 
			 np.array(data_begin_camp['offset_svrk']), YQ_i_AKNP_chan_13, 'bo', np.array(data_begin_camp['offset_svrk']), func_33, 'g')
	ax1.set_xlabel(u'offset_svrk, %')
	ax1.set_ylabel(u'$I$')
	ax1.legend((u'YQ_i_AKNP_chan_11', u'func_13', u'YQ_i_AKNP_chan_13', u'func_33'), loc='best')
	ax1.set_title(u'begin_camp.csv') 
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	savefig(r'out\csv\begin_camp.png', orientation='landscape', papertype='a4', format='png')




	






plot_1 = 0
if plot_1 == True:
	#locator = mdates.AutoDateLocator()
	#major_formatter = mdates.DateFormatter('%m-%d %H:%M')   
	plt.rcParams.update({'font.size': 14})
	if end_camp == True:
		with open(r'out\data_res_end.pkl', "rb") as resfile:
			dat = pickle.load(resfile)
	else:
		with open(r'out\data_res_begin.pkl', "rb") as resfile:
			dat = pickle.load(resfile)

	sum_k = 0
	sum_b = 0 
	parametr_W, parametr_cov = curve_fit(approk, np.array(dat['offset_svrk']), np.array(dat['YQ_W_AKNP_11']))
	#print('parametr = ', parametr)
	func_1 = parametr_W[0]*np.array(dat['offset_svrk'])+parametr_W[1]
	
	parametr_WT, parametr_cov = curve_fit(approk, np.array(dat['offset_svrk']), np.array(dat['YQ_AKPM_WT_11']))
	#print('parametr = ', parametr)
	func_2 = parametr_WT[0]*np.array(dat['offset_svrk'])+parametr_WT[1]
	sum_k += float(parametr_W[0])
	sum_k += float(parametr_WT[0])
	sum_b += float(parametr_W[1])
	sum_b += float(parametr_WT[1]) 

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(211)
	ax1.grid(True)
	ax1.plot(dat['offset_svrk'], dat['YQ_W_AKNP_11'], 'ro', dat['offset_svrk'], func_1, 'k', dat['offset_svrk'], dat['YQ_AKPM_WT_11'], 'bo', dat['offset_svrk'], func_2, 'k')
	ax1.set_xlabel(u'Офсет по СВРК, %')
	ax1.set_ylabel(u'Мощность, %$W_{ном}$')
	ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	parametr_W, parametr_cov = curve_fit(approk, np.array(dat['offset_akpm']), np.array(dat['YQ_W_AKNP_11']))
	#print('parametr = ', parametr)
	func_1 = parametr_W[0]*np.array(dat['offset_akpm'])+parametr_W[1]

	parametr_WT, parametr_cov = curve_fit(approk, np.array(dat['offset_akpm']), np.array(dat['YQ_AKPM_WT_11']))
	#print('parametr = ', parametr)
	func_2 = parametr_WT[0]*np.array(dat['offset_akpm'])+parametr_WT[1]
	sum_k += float(parametr_W[0])
	sum_k += float(parametr_WT[0])
	sum_b += float(parametr_W[1])
	sum_b += float(parametr_WT[1])

	ax2 = fig.add_subplot(212)
	ax2.grid(True)
	ax2.plot(dat['offset_akpm'], dat['YQ_W_AKNP_11'], 'ro', dat['offset_akpm'], func_1, 'k', dat['offset_akpm'], dat['YQ_AKPM_WT_11'], 'bo', dat['offset_akpm'], func_2, 'k')
	ax2.set_xlabel(u'Офсет по АКПМ, %')
	ax2.set_ylabel(u'Мощность, %$W_{ном}$')
	ax2.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#fig.autofmt_xdate()
	fig.subplots_adjust(wspace=0.05, hspace=0.15)
	if end_camp == True:
		savefig(r'out\offset_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\offset_begin.png', orientation='landscape', papertype='a4', format='png')	
	print()

	koeff_fit = sum_k/4
	koeff_fit_b = sum_b/4
	

	amperage = np.array(dat['YQ_i_AKNP_chan_11']) + (np.array(dat['YQ_i_AKNP_chan_11']) + np.array(dat['YQ_i_AKNP_chan_13'])) * koeff_fit * np.array(data_res['offset_svrk'])
	W_stok = (np.array(dat['YQ_i_AKNP_chan_11']) + np.array(dat['YQ_i_AKNP_chan_13']))/2
	W_amperage = (amperage + np.array(dat['YQ_i_AKNP_chan_13']))/2
	
	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'], W_stok, 'r', dat['time'], W_amperage, 'b')
	ax1.set_xlabel(u'Время, с')
	#ax1.set_ylabel(u'Отношение $I_{1}$ к $I_{3}$')
	ax1.legend((u'W_stok', u'W_amperage'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\W_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\W_begin.png', orientation='landscape', papertype='a4', format='png')

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(211)
	ax1.grid(True)
	ax1.plot(dat['time'], dat['YQ_W_AKNP_11'], 'r', dat['time'], dat['YQ_AKPM_WT_11'], 'b')#, dat['time'], W_stok, 'g', dat['time'], W_amperage, 'k')
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'Мощность, %$W_{ном}$')
	ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 



	ax2 = fig.add_subplot(212)
	ax2.grid(True)
	ax2.plot(dat['time'], dat['offset_svrk'], 'r', dat['time'], dat['offset_akpm'], 'b')
	ax2.set_xlabel(u'Время, с')
	ax2.set_ylabel(u'Офсет, %')
	ax2.legend((u'Офсет по СВРК,', u'Офсет по АКПМ,'), loc='best')


	#fig.autofmt_xdate()
	fig.subplots_adjust(wspace=0.05, hspace=0.15)
	if end_camp == True:
		savefig(r'out\process_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\process_begin.png', orientation='landscape', papertype='a4', format='png')	



	sum = 0 
	for i in range(len(dat['time'])):
		sum += dat['offset_svrk'][i] / dat['offset_akpm'][i]
	koeff = sum / len(dat['time'])
	
	if end_camp == True:
		print('end')
	else:
		print('begin')
	print()
	print('koeff_fit = ', koeff_fit)
	print('koeff = ', koeff)
	print()


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'], np.divide(dat['YQ_i_AKNP_chan_11'], dat['YQ_i_AKNP_chan_13']))
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'Отношение $I_{1}$ к $I_{3}$')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\amperage_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\amperage_begin.png', orientation='landscape', papertype='a4', format='png')


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'], np.array(dat['YQ_i_AKNP_chan_11']) * np.array(dat['YQ_i_AKNP_chan_13']))
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'Проивзедение $I_{1}$ и $I_{3}$')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\multiplication_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\multiplication_begin.png', orientation='landscape', papertype='a4', format='png')

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'], np.array(dat['YQ_i_AKNP_chan_11']) + np.array(dat['YQ_i_AKNP_chan_13']))
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'Сумма $I_{1}$ и $I_{3}$')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\summ_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\summ_begin.png', orientation='landscape', papertype='a4', format='png')

	
	WWW = koeff_fit*(amperage - np.array(dat['YQ_i_AKNP_chan_13']))/(amperage + np.array(dat['YQ_i_AKNP_chan_13'])) + koeff_fit_b


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'], WWW)
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'Мощность, %$W_{ном}$')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\bbb_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\bbb_begin.png', orientation='landscape', papertype='a4', format='png')
	


	dI_dOFF = []
	chan = 11
	Icp = np.sum(dat[f'YQ_i_AKNP_chan_{chan}'])/len(dat['time']) 
	for i in range(len(dat['time']) - 1): 
		dI_dOFF.append((dat[f'YQ_i_AKNP_chan_{chan}'][i+1] - dat[f'YQ_i_AKNP_chan_{chan}'][i])/(dat['offset_svrk'][i+1] - dat['offset_svrk'][i])/Icp)

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['offset_svrk'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Офсет, %')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(f'out\dI_dOFF\der\der_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\der\der_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')


	dI_dOFF = []
	chan = 12
	Icp = np.sum(dat[f'YQ_i_AKNP_chan_{chan}'])/len(dat['time']) 
	for i in range(len(dat['time']) - 1): 
		dI_dOFF.append((dat[f'YQ_i_AKNP_chan_{chan}'][i+1] - dat[f'YQ_i_AKNP_chan_{chan}'][i])/(dat['offset_svrk'][i+1] - dat['offset_svrk'][i])/Icp)

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['offset_svrk'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Офсет, %')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(f'out\dI_dOFF\der\der_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\der\der_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')


	dI_dOFF = []
	chan = 13
	Icp = np.sum(dat[f'YQ_i_AKNP_chan_{chan}'])/len(dat['time']) 
	for i in range(len(dat['time']) - 1): 
		dI_dOFF.append((dat[f'YQ_i_AKNP_chan_{chan}'][i+1] - dat[f'YQ_i_AKNP_chan_{chan}'][i])/(dat['offset_svrk'][i+1] - dat['offset_svrk'][i])/Icp)

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['offset_svrk'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Офсет, %')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\dI_dOFF_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')


	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'][1:], dI_dOFF, 'ro')
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'$dI/dOFF/Icp$, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		func ('end') 
		savefig(f'out\dI_dOFF\der\der_{chan}_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(f'out\dI_dOFF\der\der_{chan}_begin.png', orientation='landscape', papertype='a4', format='png')




	half_sum = []
	#chan = 13
	half_sum_cp = (np.sum(dat[f'YQ_i_AKNP_chan_11']) + np.sum(dat[f'YQ_i_AKNP_chan_13']))/(len(dat['time'])*2) 
	for i in range(len(dat['time']) - 1): 
		half_sum.append((((dat[f'YQ_i_AKNP_chan_11'][i+1] + dat[f'YQ_i_AKNP_chan_13'][i+1])/2) - ((dat[f'YQ_i_AKNP_chan_11'][i] + dat[f'YQ_i_AKNP_chan_13'][i])/2)) / (dat['offset_svrk'][i+1] - dat['offset_svrk'][i]) / half_sum_cp)

	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['offset_svrk'][1:], half_sum, 'ro')
	ax1.set_xlabel(u'Офсет, %')
	ax1.set_ylabel(u'dhalf_sum/dOFF/half_sum_cp, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\dhalf_sum_dOFF\dhalf_sum_dOFF_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\dhalf_sum_dOFF\dhalf_sum_dOFF_begin.png', orientation='landscape', papertype='a4', format='png')




	fig = figure(figsize=(12, 16))
	#fig.autofmt_xdate()
	ax1 = fig.add_subplot(111)
	ax1.grid(True)
	ax1.plot(dat['time'][1:], half_sum, 'ro')
	ax1.set_xlabel(u'Время, с')
	ax1.set_ylabel(u'dhalf_sum/dOFF/half_sum_cp, 1/%')
	#ax1.legend((u'Некорректированная мощность', u'Корректированная мощность'), loc='best')
	#ax1.set_xlim((dat["time"][0], dat["time"][-1]))
	#ax1.set_ylim((ylim_min, ylim_max)) 

	if end_camp == True:
		savefig(r'out\dhalf_sum_dOFF\der_end.png', orientation='landscape', papertype='a4', format='png')
	else:
		savefig(r'out\dhalf_sum_dOFF\der_begin.png', orientation='landscape', papertype='a4', format='png')














