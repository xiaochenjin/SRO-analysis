import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def plot(pSRO_jj_NN1_list,pSRO_jj_NN2_list,pSRO_jj_NN3_list,pair):
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 11}
	legend_properties = {'weight':'bold','size':11,'style': 'italic'}

	fig, axs = plt.subplots(1,1, figsize=(8,6))#,dpi=500)

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	#axs.xaxis.grid(False, which='minor')
	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))


	axs.set_xlabel("MC step",fontsize='13',weight = 'bold')
	axs.set_ylabel("SRO parameter",fontsize='13',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.0f'))

	X_axis = np.arange(1,len(pSRO_jj_NN1_list)+1)


	axs.plot(X_axis,pSRO_jj_NN1_list,color='k',linestyle='-',label='1NN {0}'.format(pair))
	axs.plot(X_axis,pSRO_jj_NN2_list,color='b',linestyle='-',label='2NN {0}'.format(pair))
	axs.plot(X_axis,pSRO_jj_NN3_list,color='g',linestyle='-',label='3NN {0}'.format(pair))

#	axs.set_xlim(-0.1,40)
#	axs.set_ylim(0,5)
	plt.legend(prop=legend_properties, frameon=False,loc='best') #bbox_to_anchor=(0.9,-0.1))#loc='best')
#   plt.title(directory_name,fontsize='13',weight = 'bold')
	plt.tight_layout()
	plt.show()

anal_dir = '/home/xcjin/Research/other-small-things/GePb/MC-sampling/128-atoms/analysis/'
N_total = 128
#c_Pb_list = ["0.0625"]
c_Pb_list = ["0.0625","0.09375","0.125"]
for conc_Pb in c_Pb_list:
	print(conc_Pb)
	c_Pb = float(conc_Pb)
	c_Ge = 1- float(conc_Pb)
	anal_conc = os.path.join(anal_dir,conc_Pb)
	infile = open(os.path.join(anal_conc,"pSRO-info-AA.txt"),'r')
	data = infile.readlines()
	energy_list = [];
	pSRO_ii_NN1_list = [];pSRO_dist_ii_NN1_list = []
	pSRO_jj_NN1_list = [];pSRO_dist_jj_NN1_list = []
	pSRO_ii_NN2_list = [];pSRO_dist_ii_NN2_list = []
	pSRO_jj_NN2_list = [];pSRO_dist_jj_NN2_list = []
	pSRO_ii_NN3_list = [];pSRO_dist_ii_NN3_list = []
	pSRO_jj_NN3_list = [];pSRO_dist_jj_NN3_list = []
	for line in data:
		energy_list.append(float(line.split()[1]))
		pSRO_ii_NN1_list.append(float(line.split()[2]))
		pSRO_ii_NN2_list.append(float(line.split()[3]))
		pSRO_ii_NN3_list.append(float(line.split()[4]))
		pSRO_jj_NN1_list.append(float(line.split()[5]))
		pSRO_jj_NN2_list.append(float(line.split()[6]))
		pSRO_jj_NN3_list.append(float(line.split()[7]))
		pSRO_dist_ii_NN1_list.append(line.split()[8])
		pSRO_dist_ii_NN2_list.append(line.split()[9])
		pSRO_dist_ii_NN3_list.append(line.split()[10])
		pSRO_dist_jj_NN1_list.append(line.split()[11])
		pSRO_dist_jj_NN2_list.append(line.split()[12])
		pSRO_dist_jj_NN3_list.append(line.split()[13])

	print("average energy of MC sampling: ", np.average(energy_list))
	print("MC sampling STD: ", np.std(energy_list))

	ave_pSRO_ii_NN1 = np.average(pSRO_ii_NN1_list)
	print("average 1NN Ge-Ge SRO parameter ", round(ave_pSRO_ii_NN1,3))
	ave_pSRO_ii_NN2 = np.average(pSRO_ii_NN2_list)
	print("average 2NN Ge-Ge SRO parameter ", round(ave_pSRO_ii_NN2,3))
	ave_pSRO_ii_NN3 = np.average(pSRO_ii_NN3_list)
	print("average 3NN Ge-Ge SRO parameter ", round(ave_pSRO_ii_NN3,3))
	ave_pSRO_jj_NN1 = np.average(pSRO_jj_NN1_list)
	print("average 1NN Pb-Pb SRO parameter ", round(ave_pSRO_jj_NN1,3))
	ave_pSRO_jj_NN2 = np.average(pSRO_jj_NN2_list)
	print("average 2NN Pb-Pb SRO parameter ", round(ave_pSRO_jj_NN2,3))
	ave_pSRO_jj_NN3 = np.average(pSRO_jj_NN3_list)
	print("average 3NN Pb-Pb SRO parameter ", round(ave_pSRO_jj_NN3,3))

	plot(pSRO_jj_NN1_list,pSRO_jj_NN2_list,pSRO_jj_NN3_list,'Pb-Pb')

	pSRO_dist_ii_NN1_list_flt = []
	for i in range(len(pSRO_dist_ii_NN1_list)):
		pSRO_dist_ii_NN1 = pSRO_dist_ii_NN1_list[i]
		pSRO_dist_ii_NN1_flt = [float(x) for x in pSRO_dist_ii_NN1.split('x')]
		pSRO_dist_ii_NN1_list_flt.append(pSRO_dist_ii_NN1_flt)

	pSRO_dist_jj_NN1_list_flt = []
	for i in range(len(pSRO_dist_jj_NN1_list)):
		pSRO_dist_jj_NN1 = pSRO_dist_jj_NN1_list[i]
		pSRO_dist_jj_NN1_flt = [float(x) for x in pSRO_dist_jj_NN1.split('x')]
		pSRO_dist_jj_NN1_list_flt.append(pSRO_dist_jj_NN1_flt)

	pSRO_dist_ii_NN2_list_flt = []
	for i in range(len(pSRO_dist_ii_NN2_list)):
		pSRO_dist_ii_NN2 = pSRO_dist_ii_NN2_list[i]
		pSRO_dist_ii_NN2_flt = [float(x) for x in pSRO_dist_ii_NN2.split('x')]
		pSRO_dist_ii_NN2_list_flt.append(pSRO_dist_ii_NN2_flt)

	pSRO_dist_jj_NN2_list_flt = []
	for i in range(len(pSRO_dist_jj_NN2_list)):
		pSRO_dist_jj_NN2 = pSRO_dist_jj_NN2_list[i]
		pSRO_dist_jj_NN2_flt = [float(x) for x in pSRO_dist_jj_NN2.split('x')]
		pSRO_dist_jj_NN2_list_flt.append(pSRO_dist_jj_NN2_flt)

	pSRO_dist_ii_NN3_list_flt = []
	for i in range(len(pSRO_dist_ii_NN3_list)):
		pSRO_dist_ii_NN3 = pSRO_dist_ii_NN3_list[i]
		pSRO_dist_ii_NN3_flt = [float(x) for x in pSRO_dist_ii_NN3.split('x')]
		pSRO_dist_ii_NN3_list_flt.append(pSRO_dist_ii_NN3_flt)

	pSRO_dist_jj_NN3_list_flt = []
	for i in range(len(pSRO_dist_jj_NN3_list)):
		pSRO_dist_jj_NN3 = pSRO_dist_jj_NN3_list[i]
		pSRO_dist_jj_NN3_flt = [float(x) for x in pSRO_dist_jj_NN3.split('x')]
		pSRO_dist_jj_NN3_list_flt.append(pSRO_dist_jj_NN3_flt)

	ave_pSRO_dist_ii_NN1 = list(np.mean(pSRO_dist_ii_NN1_list_flt,axis=0))
	ave_pSRO_dist_jj_NN1 = list(np.mean(pSRO_dist_jj_NN1_list_flt,axis=0))
	ave_pSRO_dist_ii_NN2 = list(np.mean(pSRO_dist_ii_NN2_list_flt,axis=0))
	ave_pSRO_dist_jj_NN2 = list(np.mean(pSRO_dist_jj_NN2_list_flt,axis=0))
	ave_pSRO_dist_ii_NN3 = list(np.mean(pSRO_dist_ii_NN3_list_flt,axis=0))
	ave_pSRO_dist_jj_NN3 = list(np.mean(pSRO_dist_jj_NN3_list_flt,axis=0))
	
	ave_pSRO_dist_ii_NN1 = 'x'.join([str(x) for x in ave_pSRO_dist_ii_NN1])
	ave_pSRO_dist_jj_NN1 = 'x'.join([str(x) for x in ave_pSRO_dist_jj_NN1])
	ave_pSRO_dist_ii_NN2 = 'x'.join([str(x) for x in ave_pSRO_dist_ii_NN2])
	ave_pSRO_dist_jj_NN2 = 'x'.join([str(x) for x in ave_pSRO_dist_jj_NN2])
	ave_pSRO_dist_ii_NN3 = 'x'.join([str(x) for x in ave_pSRO_dist_ii_NN3])
	ave_pSRO_dist_jj_NN3 = 'x'.join([str(x) for x in ave_pSRO_dist_jj_NN3])
		
	print("average distribution of Ge-Ge (i around i) at 1st shell: ", ave_pSRO_dist_ii_NN1)
	print("average distribution of Pb-Pb (j around j) at 1st shell: ", ave_pSRO_dist_jj_NN1)
	print("average distribution of Ge-Ge (i around i) at 2nd shell: ", ave_pSRO_dist_ii_NN2)
	print("average distribution of Pb-Pb j-i (j around j) at 2nd shell: ", ave_pSRO_dist_jj_NN2)
	print("average distribution of Ge-Ge (i around i) at 3rd shell: ", ave_pSRO_dist_ii_NN3)
	print("average distribution of Pb-Pb (j around j) at 3rd shell: ", ave_pSRO_dist_jj_NN3)

	print("\n")



