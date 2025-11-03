import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.legend_handler import HandlerTuple

def plot(file_name):
	with open (file_name) as f1:
		R_list,GR_list = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)
	
	R_kNN_list = []
	GR_kNN_list = []
	for i in range(len(R_list)):
		R = R_list[i]
		GR = GR_list[i]
		if GR > 0:
			R_kNN_list.append(R)
			GR_kNN_list.append(GR)

	N_bin = 10000
	N_atom = 64000
	a = 20
	dr = 0.5*a/N_bin
	V = a**3
	density_0 = N_atom/V
	M_list = []
	for i in range(len(GR_kNN_list)):
		R = R_kNN_list[i]
		GR = GR_kNN_list[i]
		density_R = density_0*GR
		M = density_R*(4*3.14*R*R*dr)
		M_list.append(M)

	with open ("kNN-peak-position.txt","w") as f1:
		for i in range(len(R_kNN_list)):
			f1.write(str(i+1)+'NN'+' '+str(R_kNN_list[i])+' '+str(M_list[i])+'\n')
	f1.close()

	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}

	fig, axs = plt.subplots(1,1, figsize=(6,5))#,dpi=500)

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	#axs.xaxis.grid(False, which='minor')
	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))

	axs.set_xlabel("r",fontsize='13',weight = 'bold')
	axs.set_ylabel("g(r)",fontsize='13',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

	axs.plot(R_list,GR_list,'k--')
	axs.scatter(R_kNN_list,GR_kNN_list,color='r',marker='s')
	axs.set_xlim(0,3)
#   axs.set_ylim(0,5)
	#plt.legend(prop=legend_properties, frameon=False,bbox_to_anchor=(0.9,-0.1))#loc='best')
	plt.title(file_name,fontsize='13',weight = 'bold')
	plt.show()

file_name = 'rdf-total-noperturb-before-relax-64000-atom-cell.txt'
#file_name = 'rdf-total-noperturb-before-relax-8000-atom-cell.txt'
plot(file_name)
