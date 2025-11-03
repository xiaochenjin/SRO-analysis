import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def concentrations(bond):
	if bond == "GePb" or bond == "PbPb":
		x_list = []
		for i in range(5):
			x_list.append(str(round(1 - i/(4*c_Pb),2)))
	elif bond == "GeGe" or bond == 'PbGe':
		x_list = []
		for i in range(5):
			x_list.append(str(round(1 - i/(4*c_Ge),2)))
	
	return x_list
		
		
def barplot(bond,x_list):
	#https://stackoverflow.com/questions/28750253/how-to-change-the-default-latex-font-in-matplotlib
	#https://stackoverflow.com/questions/42097053/matplotlib-cannot-find-basic-fonts
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']
	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}
	
	file_name = "{}_NN1.txt".format(bond)
	with open(file_name) as f1:
		energy_list = []
		index_list = []
		sro_list = []
		data = f1.readlines()
		for line in data:
			index = line.split()[0]
		#	if index == '0095091' or index == '0137846' or index == '0202179':
#			energy = float(line.split()[1])
			sro = line.split()[1]
#			energy_list.append(energy)
			index_list.append(index)
			sro_list.append(sro)
			sro_ave = float(line.split()[2])
			print(index)
			print("average SRO parameter: ", round(sro_ave,3))
	
	sro_list_new = []
	sro_list_int = []
	
	for i in range(len(sro_list)):
		sro_list_new.append(sro_list[i].split('x'))

	for i in range(len(sro_list_new)):
		s = sro_list_new[i]
		t = []
		for j in range(len(s)):
#			t.append(int(s[j]))
			t.append(float(s[j])/N_center) #normalized
		sro_list_int.append(t)
			
#	for i in range(len(index_list)):
#		y_list = sro_list_int[i]
#		title = str(bond)+'-NN1-'+index_list[i]

	y_list = np.mean(sro_list_int,axis = 0)
	plt.ylim(0,1)
	plt.title(title, weight = 'bold',fontsize='13')
	plt.xlabel('SRO Parameter',weight = 'bold',fontsize='13')
	plt.ylabel('Frequency',weight = 'bold',fontsize='13')
	plt.bar(x_list,y_list,label='constraints: 1NN ave',color = 'b')
	plt.legend(prop=legend_properties, frameon=False,loc='best')
	plt.show()

N_total = 128	
path = "./"
conc_Pb_list = ["0.0625","0.09375","0.125"]
bond_list = ["PbGe","GePb","PbPb"] 
title_list = ["Average Distribution (1NN Ge around Pb)","Average Distribution (1NN Pb around Ge)","Average Distribution (1NN Pb around Pb)"]

for conc in conc_Pb_list:
#	c_Ge = float(conc.split('-')[0])
#	c_Pb = float(conc.split('-')[1])
	c_Pb = conc
	c_Ge = 1 - c_Pb
	N_Ge = N_total*c_Ge
	N_Pb = N_total*c_Pb
	N_center_list = [N_Pb,N_Ge]
	for n in range(len(bond_list)):
		bond = bond_list[n]
		print(bond)
		N_center = N_center_list[n]
		x_list = concentrations(bond)
		title = title_list[n]
		barplot(bond,x_list)
