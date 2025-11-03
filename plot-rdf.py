import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.legend_handler import HandlerTuple
from scipy.signal import savgol_filter

def smooth_GR(GR_list):
	window_length = 5
	smoothed_GR_list = savgol_filter(GR_list,window_length,polyorder=1)
	return smoothed_GR_list

def plot(file_name):
	with open (file_name) as f1:
		R_list,GR_list = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)
	smoothed_GR_list = smooth_GR(GR_list)

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

	axs.set_xlabel("r (angstrom)",fontsize='13',weight = 'bold')
	axs.set_ylabel("g(r)",fontsize='13',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

	axs.plot(R_list,[1]*len(R_list),'b--')
	axs.scatter(R_list,GR_list,color='r',marker='o',s=1,alpha = 1,label='raw APT data')
#	axs.plot(R_list,smoothed_GR_list,'r-',label = 'smoothed RDF')
	axs.set_xlim(-0.1,5)
	axs.set_ylim(0,5)
	plt.legend(prop=legend_properties, frameon=False,loc='best')
#	plt.title(file_name,fontsize='13',weight = 'bold')
	plt.tight_layout()
	plt.show()

cell_dimension = int(input("Dimension (nm): "))
file_name = "rdf-{0}x{0}x{0}-cube-GeSn.txt".format(str(cell_dimension))
#file_name = 'rdf-specific-cube-GeSn.txt'
plot(file_name)
