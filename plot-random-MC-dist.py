import os
import numpy as np
import matplotlib.pyplot as plt

N_total = 128
anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'
bond_type = 'AA'

# concentration of neighboring atoms
c_Pb_list = ["0.0625","0.09375","0.125"]

def barplot(x_list,y_MC_list,y_RM_list,title):
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']
	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}

	ind = np.arange(len(x_list))
	width = 0.25

	plt.ylim(0,1)
	plt.title(title,weight = 'bold',fontsize='11')
#   plt.title('theoretical distribution for random alloy',weight = 'bold',fontsize='13')
	plt.xlabel('SRO Parameter',weight = 'bold',fontsize='11')
	plt.ylabel('Frequency',weight = 'bold',fontsize='11')
	#plt.xticks(fontproperties)
	#plt.yticks(fontproperties)
	#plt.bar(x_list,y_list2,label=label_list[i],color = 'r')
	plt.bar(ind,y_MC_list,width,label= 'average of MC sampling',color = 'black')
	plt.bar(ind+width,y_RM_list,width,label= 'theoretical distribution of random alloy',color = 'red')
	plt.xticks(ind+width,x_list,fontsize ='7',weight = 'bold')
	plt.legend(prop=legend_properties, frameon=False,loc='best')
	plt.show()

for i in range(len(c_Pb_list)):
	print(c_Pb_list[i])
	anal_conc = os.path.join(anal_dir,c_Pb_list[i])
	conc_j = float(c_Pb_list[i])
	conc_i = 1 - conc_j

	x_ii_1NN_list = [] #i around i
	for j in range(5):
		x_ii_1NN_list.append(str(round(1 - j/(4*conc_i),2)))
	
	x_jj_1NN_list = [] #j around j
	for j in range(5):
		x_jj_1NN_list.append(str(round(1 - j/(4*conc_j),2)))

	x_ii_2NN_list = [] #i around i
	for j in range(13):
		x_ii_2NN_list.append(str(round(1 - j/(12*conc_i),2)))

	x_jj_2NN_list = [] #j around j
	for j in range(13):
		x_jj_2NN_list.append(str(round(1 - j/(12*conc_j),2)))


	x_ii_3NN_list = x_ii_2NN_list
	x_jj_3NN_list = x_jj_2NN_list

	x_ij_1NN_list = x_jj_1NN_list #j around i
	x_ji_1NN_list = x_ii_1NN_list #i around j
	x_ij_2NN_list = x_jj_2NN_list #j around i
	x_ji_2NN_list = x_ii_2NN_list #i around j
	x_ij_3NN_list = x_jj_3NN_list #j around i
	x_ji_3NN_list = x_ii_3NN_list #i around j

	infile =  open(os.path.join(anal_conc,"pSRO-ave-MC-{}.txt".format(bond_type)),'r')
	data = infile.readlines()
	infile.close()
	ave_jj_MC_1NN = float(data[5].split()[-1])
	ave_jj_MC_2NN = float(data[6].split()[-1])
	dist_jj_MC_1NN = list(map(float,data[9].split()[-1].split('x')))
	dist_jj_MC_2NN = list(map(float,data[11].split()[-1].split('x')))
	N_center_jj = N_total*conc_j
	dist_jj_MC_1NN = [x/N_center_jj for x in dist_jj_MC_1NN]
	dist_jj_MC_2NN = [x/N_center_jj for x in dist_jj_MC_2NN]
	print("average 1NN Pb-Pb SRO parameter (MC): ", ave_jj_MC_1NN)
	print("average 2NN Pb-Pb SRO parameter (MC): ", ave_jj_MC_2NN)

	infile =  open(os.path.join(anal_conc,"pSRO-theoretical-random-{}.txt".format(bond_type)),'r')
	data = infile.readlines()
	infile.close()
	dist_jj_RM_1NN = list(map(float,data[7].split()[-1].split('x')))
	dist_jj_RM_2NN = list(map(float,data[9].split()[-1].split('x')))
	N_center_jj = N_total*conc_j

	#double check correction of x-axis:"
	ave_jj_1NN_SRO_MC = 0
	for j in range(5):
		ave_jj_1NN_SRO_MC += float(x_jj_1NN_list[j])*dist_jj_RM_1NN[j]
	print("average 1NN Pb-Pb SRO parameter (random): ", ave_jj_1NN_SRO_MC)

	ave_jj_2NN_SRO_MC = 0
	for j in range(13):
		ave_jj_2NN_SRO_MC += float(x_jj_2NN_list[j])*dist_jj_RM_2NN[j]
	print("average 2NN Pb-Pb SRO parameter (random): ", ave_jj_2NN_SRO_MC)

	barplot(x_jj_1NN_list,dist_jj_MC_1NN,dist_jj_RM_1NN,"Pb-{}: 1NN Pb-Pb SRO parameters".format(conc_j))
	barplot(x_jj_2NN_list,dist_jj_MC_2NN,dist_jj_RM_2NN,"Pb-{}: 2NN Pb-Pb SRO parameters".format(conc_j))
		


