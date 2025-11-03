import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

N_total = 128
anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'
bond_type = 'AA'

# concentration of neighboring atoms
def ave_plot(conc_list,ave_SRO_list,title):
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']
	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}

	plt.axes().yaxis.set_minor_locator(AutoMinorLocator(2))
	plt.axes().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.ylim(-1,1)
	plt.title(title,weight = 'bold',fontsize='13')
	#   plt.title('theoretical distribution for random alloy',weight = 'bold',fontsize='13')
	plt.xlabel('Pb concentration',weight = 'bold',fontsize='13')
	plt.ylabel('SRO parameter',weight = 'bold',fontsize='13')
	plt.xticks(weight='bold',fontsize='11')
	plt.yticks(weight='bold',fontsize='11')
	plt.plot(conc_list,ave_SRO_list,'ro-',label= 'average of MC sampling')
	line_RM = [0]*len(conc_list)
	plt.plot(conc_list,line_RM, 'k--',label= 'random alloy')
	plt.legend(prop=legend_properties, frameon=False,loc='best')
	plt.tight_layout()
	plt.show()

c_Pb_list = ["0.0625","0.09375","0.125"]
ave_jj_MC_1NN_list = []; ave_jj_MC_2NN_list = []
conc_j_list = []
for i in range(len(c_Pb_list)):
	print(c_Pb_list[i])
	anal_conc = os.path.join(anal_dir,c_Pb_list[i])
	conc_j = float(c_Pb_list[i])
	conc_i = 1 - conc_j
	conc_j_list.append(conc_j)

	infile =  open(os.path.join(anal_conc,"pSRO-ave-MC-{}.txt".format(bond_type)),'r')
	data = infile.readlines()
	infile.close()
	ave_jj_MC_1NN = float(data[5].split()[-1])
	ave_jj_MC_2NN = float(data[6].split()[-1])

	ave_jj_MC_1NN_list.append(ave_jj_MC_1NN)
	ave_jj_MC_2NN_list.append(ave_jj_MC_2NN)

ave_plot(conc_j_list,ave_jj_MC_1NN_list,"Average 1NN Pb-Pb SRO Parameter")
ave_plot(conc_j_list,ave_jj_MC_2NN_list,"Average 2NN Pb-Pb SRO Parameter")		

