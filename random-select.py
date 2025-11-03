import os
import numpy as np
import shutil
import random
import operator

anal_dir = '/home/xcjin/Research/other-small-things/GePb/MC-sampling/128-atoms/analysis/'
data_dir = '/home/xcjin/Research/other-small-things/GePb/MC-sampling/128-atoms/data'

equil = input("Equilibrium or not? ")
N_random = int(input("Number of random configurations to select for each concentration: "))

c_Pb_list = ["0.0625","0.09375","0.125"]
for c_Pb in c_Pb_list:
	anal_conc = os.path.join(anal_dir,c_Pb)
	data_conc = os.path.join(data_dir,c_Pb)

	infile = open(os.path.join(anal_conc,'pSRO-info-AA.txt'),'r')
	if equil == 'N': data=infile.readlines()[:100]
	if equil == 'Y': data=infile.readlines()[300:]
	infile.close()
	index_list = []
	energy_list = [];
	pSRO_jj_NN1_list = [];pSRO_dist_jj_NN1_list = []
	pSRO_jj_NN2_list = [];pSRO_dist_jj_NN2_list = []
	pSRO_jj_NN3_list = [];pSRO_dist_jj_NN3_list = []
	for line in data:
		if not line.split()[0] in index_list:
			index_list.append(line.split()[0])
			energy_list.append(float(line.split()[1]))
			pSRO_jj_NN1_list.append(float(line.split()[5]))
			pSRO_jj_NN2_list.append(float(line.split()[6]))
			pSRO_jj_NN3_list.append(float(line.split()[7]))
			pSRO_dist_jj_NN1_list.append(line.split()[11])
			pSRO_dist_jj_NN2_list.append(line.split()[12])
			pSRO_dist_jj_NN3_list.append(line.split()[13])
	
	infor_totlist = list(zip(index_list,energy_list,pSRO_jj_NN1_list,pSRO_jj_NN2_list,pSRO_jj_NN3_list,pSRO_dist_jj_NN1_list,pSRO_dist_jj_NN2_list,pSRO_dist_jj_NN3_list))

	
	sampled_infor_list = random.sample(infor_totlist,N_random)
	sampled_infor_list.sort(key=operator.itemgetter(0))
	
	copydir = os.path.join(anal_conc,"randomly-select")
	if not os.path.exists(copydir):
		os.mkdir(copydir)

	if equil == 'N': 
		folder_name = "before-equilibrium"
		file_name = "randomly-selected-structures-before-equilibrium.txt"

	if equil == 'Y':
		folder_name = "after-equilibrium"
		file_name = "randomly-selected-structures-after-equilibrium.txt"

	if os.path.exists(os.path.join(copydir,folder_name)):
		shutil.rmtree(os.path.join(copydir,folder_name))
	os.mkdir(os.path.join(copydir,folder_name))

	with open (os.path.join(copydir,file_name),'w') as f1:
		for i in range(len(sampled_infor_list)):
			infor = sampled_infor_list[i]
			index = infor[0]
			shutil.copyfile(os.path.join(data_conc,"contcar-only/CONTCAR-{0}".format(index)),os.path.join(os.path.join(copydir,folder_name),"CONTCAR-{0}".format(index)))				
			f1.write(str(infor[0])+' '+str(infor[1])+' '+str(infor[2])+' '
					+str(infor[3])+' '+str(infor[4])+' '+str(infor[5])+' '
					+str(infor[6])+' '+str(infor[7])+'\n')
	f1.close()


