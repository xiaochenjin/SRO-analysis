import numpy as np
import os
import computeSRO_dist

N_species = 2
N_shell = 3
num_cell = [4,4,4]
N_cell = num_cell[0]*num_cell[1]*num_cell[2]
type_pair = 'A-A'

#Ge-Pb:boundaries for AB pairs
conc_Pb_list = ['0.0625','0.09375','0.125']
BC_GePb_CN_list=[[2.59,2.67],[2.59,2.69],[2.61,2.70]]
BC_GePb_NN2_list=[[3.88,4.20],[3.81,4.25],[3.83,4.31]]
BC_GePb_NN3_list=[[4.26,4.92],[4.26,5.06],[4.31,5.23]]

BC_AB_CN_totlist = [BC_GePb_CN_list]
BC_AB_NN2_totlist = [BC_GePb_NN2_list]
BC_AB_NN3_totlist = [BC_GePb_NN3_list]

BC_AB_totlist = [BC_AB_CN_totlist,BC_AB_NN2_totlist,BC_AB_NN3_totlist]

#Ge-Ge and Pb-Pb boundaries (AA pairs)
BC_PbPb_CN_list=[[2.77,2.81],[2.78,2.82],[2.78,2.87]]
BC_PbPb_NN2_list=[[4.05,4.16],[4.05,4.20],[3.9,4.27]]
BC_PbPb_NN3_list=[[4.55,4.79],[4.57,4.89],[4.36,4.95]]

BC_GeGe_CN_list=[[2.37,2.50],[2.36,2.51],[2.38,2.54]]
BC_GeGe_NN2_list=[[3.72,4.49],[3.72,4.52],[3.70,4.54]]
BC_GeGe_NN3_list=[[4.49,5.12],[4.52,5.17],[4.54,5.26]]

BC_AA_CN_totlist = [BC_GeGe_CN_list,BC_PbPb_CN_list]
BC_AA_NN2_totlist = [BC_GeGe_NN2_list,BC_PbPb_NN2_list]
BC_AA_NN3_totlist = [BC_GeGe_NN3_list,BC_PbPb_NN3_list]

BC_AA_totlist = [BC_AA_CN_totlist,BC_AA_NN2_totlist,BC_AA_NN3_totlist]

anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'
data_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/data/'
for conc_Pb in conc_Pb_list:
	print(conc_Pb)
	data_conc = os.path.join(data_dir,conc_Pb)
	anal_conc = os.path.join(anal_dir,conc_Pb)
	infile = open(os.path.join(data_conc,'process_doublecount.txt'),'r')
	data = infile.readlines()
	index_list = [];energy_list = []
	for line in data:
		if "Accepted" in line:
			index = line.split()[0]
			energy = line.split()[1]
			index_list.append(index)
			energy_list.append(energy)

	pSRO_actual_list = []
	pSRO_dist_totlist = []
	for index in index_list:
		pSRO_actual,pSRO_std_actual,pSRO_dist_actual = computeSRO_dist.PSRO(index,N_species,N_shell,N_cell,type_pair,data_dir,conc_Pb,conc_Pb_list,BC_AA_totlist,BC_AB_totlist)
		pSRO_dist_actual_write_list = []
		for i in range(len(pSRO_dist_actual)):
			pSRO_dist_actual_sgl = pSRO_dist_actual[i] #SRO for one A-A pair (for SiGeSn, Si-Si, Ge-Ge, Sn-Sn))
			for j in range(len(pSRO_dist_actual_sgl)):
				pSRO_dist_actual_shell = pSRO_dist_actual_sgl[j] #SRO for one shell (i-j and j-i, 1NN, 2NN...)
				pSRO_dist_actual_write_shell = 'x'.join(map(str,pSRO_dist_actual_shell))
				pSRO_dist_actual_write_list.append(pSRO_dist_actual_write_shell)
		pSRO_actual_list.append(pSRO_actual)
		pSRO_dist_totlist.append(pSRO_dist_actual_write_list)
	#	print(index)
	#	print(pSRO_actual)
	#	print(pSRO_dist_actual_write_list)

	with open(os.path.join(anal_conc,'pSRO-info-AA.txt'),'w') as f1:
		for i in range(len(index_list)):
			f1.write(str(index_list[i])+' '+str(energy_list[i])+' ')
			pSRO_actual = pSRO_actual_list[i]
			for j in range(len(pSRO_actual)):
				f1.write(str(pSRO_actual[j])+' ')
			pSRO_dist = pSRO_dist_totlist[i]
			for j in range(len(pSRO_dist)):
				f1.write(str(pSRO_dist[j])+' ')
			f1.write('\n')
	f1.close()

