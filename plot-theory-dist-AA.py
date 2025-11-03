import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb

N_total = 128
anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'

# concentration of neighboring atoms
c_Pb_list = ["0.0625","0.09375","0.125"]

def barplot(x_list,y_list,title):
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']
	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}

	plt.ylim(0,1)
	plt.title(title,weight = 'bold',fontsize='13')
#   plt.title('theoretical distribution for random alloy',weight = 'bold',fontsize='13')
	plt.xlabel('SRO Parameter',weight = 'bold',fontsize='13')
	plt.ylabel('Frequency',weight = 'bold',fontsize='13')
	#plt.xticks(fontproperties)
	#plt.yticks(fontproperties)
	#plt.bar(x_list,y_list2,label=label_list[i],color = 'r')
	plt.bar(x_list,y_list,label='random alloy',color = 'r')
	plt.legend(prop=legend_properties, frameon=False,loc='best')
	plt.show()

for i in range(len(c_Pb_list)):
	print(c_Pb_list[i])
	anal_conc = os.path.join(anal_dir,c_Pb_list[i])
	conc_j = float(c_Pb_list[i])
	conc_i = 1 - conc_j
	#1st NN;
	title_list = ['Theoretical Distribution (1NN Ge around Ge)','Theoretical Distribution (1NN Pb around Pb)']
	
	x_ii_1NN_list = [] #i around i
	for j in range(5):
		x_ii_1NN_list.append(str(round(1 - j/(4*conc_i),2)))
	x_ii_1NN_list_flt = []
	for j in range(5):
		x_ii_1NN_list_flt.append(1 - j/(4*conc_i))
	y1 = conc_j**4 #occupacy  = 0
	y2 = 4*(conc_j**3)*conc_i #occupacy  = 1 
	y3 = 6*(conc_j**2)*(conc_i**2) #occupacy  = 2
	y4 = 4*conc_j*conc_i**3 #occupacy  = 3
	y5 = conc_i**4 # occupacy = 4
	y_ii_1NN_list = [y1,y2,y3,y4,y5]
#	print(y_list)
	N_center_ii = N_total*conc_i #center atom: the frequency counts number of centered atom that have SRO parameter of ...
	y_ii_1NN_list2 = [N_center_ii * y for y in y_ii_1NN_list]
	ave_ii_1NN_SRO = 0
	for j in range(5):
		ave_ii_1NN_SRO += x_ii_1NN_list_flt[j]*y_ii_1NN_list[j]
	print("average SRO parameter (1NN i around i): ", ave_ii_1NN_SRO)

	x_jj_1NN_list = [] #j around j
	for j in range(5):
		x_jj_1NN_list.append(str(round(1 - j/(4*conc_j),2)))
	x_jj_1NN_list_flt = []
	for j in range(5):
		x_jj_1NN_list_flt.append(1 - j/(4*conc_j))
	
	y_jj_1NN_list = np.flip(y_ii_1NN_list)
	N_center_jj = N_total*conc_j
	y_jj_1NN_list2 = [N_center_jj * y for y in y_jj_1NN_list]
	ave_jj_1NN_SRO = 0
	for j in range(5):
		ave_jj_1NN_SRO += x_jj_1NN_list_flt[j]*y_jj_1NN_list[j]
	print("average SRO parameter (1NN j around j): ", ave_jj_1NN_SRO)

	#barplot(x_ii_1NN_list,y_ii_1NN_list,title_list[0])
	#barplot(x_jj_1NN_list,y_jj_1NN_list,title_list[1])

	#2nd NN;
	title_list = ['Theoretical Distribution (2NN Ge around Ge)','Theoretical Distribution (2NN Pb around Pb)']

	x_ii_2NN_list = [] #i around i
	for j in range(13):
		x_ii_2NN_list.append(str(round(1 - j/(12*conc_i),2)))
	x_ii_2NN_list_flt = []
	for j in range(13):
		x_ii_2NN_list_flt.append(1 - j/(12*conc_i))
	y_ii_2NN_list = []
	for occupancy in range(13):
		y = comb(12,occupancy)*(conc_i**occupancy)*(conc_j**(12-occupancy))
		y_ii_2NN_list.append(y)
#   print(y_list)
	N_center_ii = N_total*conc_i #center atom: the frequency counts number of centered atom that have SRO parameter of ...
	y_ii_2NN_list2 = [N_center_ii * y for y in y_ii_2NN_list]
	ave_ii_2NN_SRO = 0
	for j in range(13):
		ave_ii_2NN_SRO += x_ii_2NN_list_flt[j]*y_ii_2NN_list[j]
	print("average SRO parameter (2NN i around i): ", ave_ii_2NN_SRO)


	x_jj_2NN_list = [] #j around j
	for j in range(13):
		x_jj_2NN_list.append(str(round(1 - j/(12*conc_j),2)))
	x_jj_2NN_list_flt = []
	for j in range(13):
		x_jj_2NN_list_flt.append(1 - j/(12*conc_j))
	y_jj_2NN_list = np.flip(y_ii_2NN_list)
	N_center_jj = N_total*conc_j
	y_jj_2NN_list2 = [N_center_jj * y for y in y_jj_2NN_list]
	ave_jj_2NN_SRO = 0
	for j in range(13):
		ave_jj_2NN_SRO += x_jj_2NN_list_flt[j]*y_jj_2NN_list[j]
	print("average SRO parameter (2NN j around j): ", ave_jj_2NN_SRO)

	#barplot(x_ii_2NN_list,y_ii_2NN_list,title_list[0])
	#barplot(x_jj_2NN_list,y_jj_2NN_list,title_list[1])

	#3rd NN;
	title_list = ['Theoretical Distribution (3NN Ge around Ge)','Theoretical Distribution (3NN Pb around Pb)']

	x_ii_3NN_list = [] #i around i
	for j in range(13):
		x_ii_3NN_list.append(str(round(1 - j/(12*conc_i),2)))
	x_ii_3NN_list_flt = []
	for j in range(13):
		x_ii_3NN_list_flt.append(1 - j/(12*conc_i))
	y_ii_3NN_list = []
	for occupancy in range(13):
		y = comb(12,occupancy)*(conc_i**occupancy)*(conc_j**(12-occupancy))
		y_ii_3NN_list.append(y)
#   print(y_list)
	N_center_ii = N_total*conc_i #center atom: the frequency counts number of centered atom that have SRO parameter of ...
	y_ii_3NN_list2 = [N_center_ii * y for y in y_ii_3NN_list]
	ave_ii_3NN_SRO = 0
	for j in range(13):
		ave_ii_3NN_SRO += x_ii_3NN_list_flt[j]*y_ii_3NN_list[j]
	print("average SRO parameter (3NN i around i): ", ave_ii_3NN_SRO)

	x_jj_3NN_list = [] #j around j
	for j in range(13):
		x_jj_3NN_list.append(str(round(1 - j/(12*conc_j),2)))
	x_jj_3NN_list_flt = []
	for j in range(13):
		x_jj_3NN_list_flt.append(1 - j/(12*conc_j))
	y_jj_3NN_list = np.flip(y_ii_3NN_list)
	N_center_jj = N_total*conc_j
	y_jj_3NN_list2 = [N_center_jj * y for y in y_jj_3NN_list]
	ave_jj_3NN_SRO = 0
	for j in range(13):
		ave_jj_3NN_SRO += x_jj_3NN_list_flt[j]*y_jj_3NN_list[j]
	print("average SRO parameter (3NN j around j): ", ave_jj_3NN_SRO)

	#barplot(x_ii_3NN_list,y_ii_3NN_list,title_list[0])
	#barplot(x_jj_3NN_list,y_jj_3NN_list,title_list[1])

	y_ii_1NN_list = 'x'.join([str(x) for x in y_ii_1NN_list])
	y_jj_1NN_list = 'x'.join([str(x) for x in y_jj_1NN_list])
	y_ii_2NN_list = 'x'.join([str(x) for x in y_ii_2NN_list])
	y_jj_2NN_list = 'x'.join([str(x) for x in y_jj_2NN_list])
	y_ii_3NN_list = 'x'.join([str(x) for x in y_ii_3NN_list])
	y_jj_3NN_list = 'x'.join([str(x) for x in y_jj_3NN_list])


	with open(os.path.join(anal_conc,"pSRO-theoretical-random-AA.txt"),"w") as f1:
		f1.write("average 1NN Ge-Ge SRO parameter: "+str(round(ave_ii_1NN_SRO,3))+'\n')
		f1.write("average 2NN Ge-Ge SRO parameter: "+str(round(ave_ii_2NN_SRO,3))+'\n')
		f1.write("average 3NN Ge-Ge SRO parameter: "+str(round(ave_ii_3NN_SRO,3))+'\n')
		f1.write("average 1NN Pb-Pb SRO parameter: "+str(round(ave_jj_1NN_SRO,3))+'\n')
		f1.write("average 2NN Pb-Pb SRO parameter: "+str(round(ave_jj_2NN_SRO,3))+'\n')
		f1.write("average 3NN Pb-Pb SRO parameter: "+str(round(ave_jj_3NN_SRO,3))+'\n')
		f1.write("average distribution of Ge-Ge (i around i) at 1st shell: "+str(y_ii_1NN_list)+'\n')
		f1.write("average distribution of Pb-Pb (j around j) at 1st shell: "+str(y_jj_1NN_list)+'\n')
		f1.write("average distribution of Ge-Ge (i around i) at 2nd shell: "+str(y_ii_2NN_list)+'\n')
		f1.write("average distribution of Pb-Pb (j around j) at 2nd shell: "+str(y_jj_2NN_list)+'\n')
		f1.write("average distribution of Ge-Ge (i around i) at 3rd shell: "+str(y_ii_3NN_list)+'\n')
		f1.write("average distribution of Pb-Pb (j around j) at 3rd shell: "+str(y_jj_3NN_list)+'\n')
	f1.close()	
