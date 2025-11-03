import os
import numpy as np
import matplotlib.pyplot as plt

def concentrations(bond):
	if bond == "SiSn" or bond == "SnSn":
		x_list = []
		for i in range(13):
			x_list.append(str(round(1 - i/(4*c_Sn),2)))
	elif bond == "SiSi":
		x_list = []
		for i in range(13):
			x_list.append(str(round(1 - i/(4*c_Si),2)))
	
	return x_list
		
		
def barplot(bond,x_list):
	#https://stackoverflow.com/questions/28750253/how-to-change-the-default-latex-font-in-matplotlib
	#https://stackoverflow.com/questions/42097053/matplotlib-cannot-find-basic-fonts
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']
	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}
	
	file_name = "{}_NN3.txt".format(bond)
	with open(file_name) as f1:
		energy_list = []
		index_list = []
		sro_list = []
		data = f1.readlines()
		for line in data:
			index = line.split()[0]
			if index == '0095091' or index == '0137846' or index == '0202179':
				energy = float(line.split()[1])
				sro = line.split()[2]
				energy_list.append(energy)
				index_list.append(index)
				sro_list.append(sro)
	
	sro_list_new = []
	sro_list_int = []
	
	for i in range(len(sro_list)):
		sro_list_new.append(sro_list[i].split('x'))

	for i in range(len(sro_list_new)):
		s = sro_list_new[i]
		t = []
		for j in range(len(s)):
			t.append(int(s[j]))
		sro_list_int.append(t)
			
	for i in range(len(index_list)):
		y_list = sro_list_int[i]
		title = str(bond)+'-NN3-'+index_list[i]
		plt.title(title)
		plt.xlabel('SRO Parameter',weight = '13',fontsize='13')
		plt.ylabel('Frequency',weight = '13',fontsize='13')
		#plt.xticks(fontproperties)
		#plt.yticks(fontproperties)
		plt.bar(x_list,y_list)
		plt.show()
	
path = "./"
conc_list = ["0.75-0.25"]
bond_list = ["SiSi","SiSn","SnSn"]

for conc in conc_list:
	c_Si = float(conc.split('-')[0])
	c_Sn = float(conc.split('-')[1])
	
	for bond in bond_list:
		x_list = concentrations(bond)
		barplot(bond,x_list)
