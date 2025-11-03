import os
import numpy as np

anal_dir = '/home/xcjin/Research/other-small-things/GePb/MC-sampling/128-atoms/analysis/'
N_total = 128
N_exclude = 100
c_Pb_list = ["0.0625","0.09375","0.125"]
for conc_Pb in c_Pb_list:
	print(conc_Pb)
	c_Pb = float(conc_Pb)
	c_Ge = 1- float(conc_Pb)
	anal_conc = os.path.join(anal_dir,conc_Pb)
	infile = open(os.path.join(anal_conc,"pSRO-info-AB.txt"),'r')
	data = infile.readlines()[N_exclude:]
	energy_list = [];
	pSRO_NN1_list = [];pSRO_dist_ij_NN1_list = [];pSRO_dist_ji_NN1_list = []
	pSRO_NN2_list = [];pSRO_dist_ij_NN2_list = [];pSRO_dist_ji_NN2_list = []
	pSRO_NN3_list = [];pSRO_dist_ij_NN3_list = [];pSRO_dist_ji_NN3_list = []
	for line in data:
		energy_list.append(float(line.split()[1]))
		pSRO_NN1_list.append(float(line.split()[2]))
		pSRO_NN2_list.append(float(line.split()[3]))
		pSRO_NN3_list.append(float(line.split()[4]))
		pSRO_dist_ij_NN1_list.append(line.split()[5])
		pSRO_dist_ji_NN1_list.append(line.split()[6])
		pSRO_dist_ij_NN2_list.append(line.split()[7])
		pSRO_dist_ji_NN2_list.append(line.split()[8])
		pSRO_dist_ij_NN3_list.append(line.split()[9])
		pSRO_dist_ji_NN3_list.append(line.split()[10])

	print("average energy of MC sampling: ", np.average(energy_list))
	print("MC sampling STD: ", np.std(energy_list))

	ave_pSRO_NN1 = np.average(pSRO_NN1_list)
	ave_pSRO_NN2 = np.average(pSRO_NN2_list)
	ave_pSRO_NN3 = np.average(pSRO_NN3_list)
	print("average 1NN Ge-Pb SRO parameter ", round(ave_pSRO_NN1,3))
	print("average 2NN Ge-Pb SRO parameter ", round(ave_pSRO_NN2,3))
	print("average 3NN Ge-Pb SRO parameter ", round(ave_pSRO_NN3,3))

	pSRO_dist_ij_NN1_list_flt = []
	for i in range(len(pSRO_dist_ij_NN1_list)):
		pSRO_dist_ij_NN1 = pSRO_dist_ij_NN1_list[i]
		pSRO_dist_ij_NN1_flt = [float(x) for x in pSRO_dist_ij_NN1.split('x')]
		pSRO_dist_ij_NN1_list_flt.append(pSRO_dist_ij_NN1_flt)

	pSRO_dist_ji_NN1_list_flt = []
	for i in range(len(pSRO_dist_ji_NN1_list)):
		pSRO_dist_ji_NN1 = pSRO_dist_ji_NN1_list[i]
		pSRO_dist_ji_NN1_flt = [float(x) for x in pSRO_dist_ji_NN1.split('x')]
		pSRO_dist_ji_NN1_list_flt.append(pSRO_dist_ji_NN1_flt)

	pSRO_dist_ij_NN2_list_flt = []
	for i in range(len(pSRO_dist_ij_NN2_list)):
		pSRO_dist_ij_NN2 = pSRO_dist_ij_NN2_list[i]
		pSRO_dist_ij_NN2_flt = [float(x) for x in pSRO_dist_ij_NN2.split('x')]
		pSRO_dist_ij_NN2_list_flt.append(pSRO_dist_ij_NN2_flt)

	pSRO_dist_ji_NN2_list_flt = []
	for i in range(len(pSRO_dist_ji_NN2_list)):
		pSRO_dist_ji_NN2 = pSRO_dist_ji_NN2_list[i]
		pSRO_dist_ji_NN2_flt = [float(x) for x in pSRO_dist_ji_NN2.split('x')]
		pSRO_dist_ji_NN2_list_flt.append(pSRO_dist_ji_NN2_flt)

	pSRO_dist_ij_NN3_list_flt = []
	for i in range(len(pSRO_dist_ij_NN3_list)):
		pSRO_dist_ij_NN3 = pSRO_dist_ij_NN3_list[i]
		pSRO_dist_ij_NN3_flt = [float(x) for x in pSRO_dist_ij_NN3.split('x')]
		pSRO_dist_ij_NN3_list_flt.append(pSRO_dist_ij_NN3_flt)

	pSRO_dist_ji_NN3_list_flt = []
	for i in range(len(pSRO_dist_ji_NN3_list)):
		pSRO_dist_ji_NN3 = pSRO_dist_ji_NN3_list[i]
		pSRO_dist_ji_NN3_flt = [float(x) for x in pSRO_dist_ji_NN3.split('x')]
		pSRO_dist_ji_NN3_list_flt.append(pSRO_dist_ji_NN3_flt)

	ave_pSRO_dist_ij_NN1 = list(np.mean(pSRO_dist_ij_NN1_list_flt,axis=0))
	ave_pSRO_dist_ji_NN1 = list(np.mean(pSRO_dist_ji_NN1_list_flt,axis=0))
	ave_pSRO_dist_ij_NN2 = list(np.mean(pSRO_dist_ij_NN2_list_flt,axis=0))
	ave_pSRO_dist_ji_NN2 = list(np.mean(pSRO_dist_ji_NN2_list_flt,axis=0))
	ave_pSRO_dist_ij_NN3 = list(np.mean(pSRO_dist_ij_NN3_list_flt,axis=0))
	ave_pSRO_dist_ji_NN3 = list(np.mean(pSRO_dist_ji_NN3_list_flt,axis=0))

	ave_pSRO_dist_ij_NN1 = 'x'.join([str(x) for x in ave_pSRO_dist_ij_NN1])
	ave_pSRO_dist_ji_NN1 = 'x'.join([str(x) for x in ave_pSRO_dist_ji_NN1])
	ave_pSRO_dist_ij_NN2 = 'x'.join([str(x) for x in ave_pSRO_dist_ij_NN2])
	ave_pSRO_dist_ji_NN2 = 'x'.join([str(x) for x in ave_pSRO_dist_ji_NN2])
	ave_pSRO_dist_ij_NN3 = 'x'.join([str(x) for x in ave_pSRO_dist_ij_NN3])
	ave_pSRO_dist_ji_NN3 = 'x'.join([str(x) for x in ave_pSRO_dist_ji_NN3])
		
	print("average distribution of i-j (j around i) at 1st shell: ", ave_pSRO_dist_ij_NN1)
	print("average distribution of j-i (i around j) at 1st shell: ", ave_pSRO_dist_ji_NN1)
	print("average distribution of i-j (j around i) at 2nd shell: ", ave_pSRO_dist_ij_NN2)
	print("average distribution of j-i (i around j) at 2nd shell: ", ave_pSRO_dist_ji_NN2)
	print("average distribution of i-j (j around i) at 3rd shell: ", ave_pSRO_dist_ij_NN3)
	print("average distribution of j-i (i around j) at 3rd shell: ", ave_pSRO_dist_ji_NN3)

	print("\n")

	with open (os.path.join(anal_conc,"pSRO-ave-MC-AB.txt"),"w") as f1:
		f1.write("average energy of MC sampling: "+str(np.average(energy_list))+'\n')
		f1.write("MC sampling STD: "+str(np.std(energy_list))+'\n')
		f1.write("average 1NN Ge-Pb SRO parameter: "+str(round(ave_pSRO_NN1,3))+'\n')
		f1.write("average 2NN Ge-Pb SRO parameter: "+str(round(ave_pSRO_NN2,3))+'\n')
		f1.write("average 3NN Ge-Pb SRO parameter: "+str(round(ave_pSRO_NN3,3))+'\n')
		f1.write("average distribution of Ge-Pb (j around i) at 1st shell: "+str(ave_pSRO_dist_ij_NN1)+'\n')
		f1.write("average distribution of Pb-Ge (i around j) at 1st shell: "+str(ave_pSRO_dist_ji_NN1)+'\n')
		f1.write("average distribution of Ge-Pb (j around i) at 2nd shell: "+str(ave_pSRO_dist_ij_NN2)+'\n')
		f1.write("average distribution of Pb-Ge (i around j) at 2nd shell: "+str(ave_pSRO_dist_ji_NN2)+'\n')
		f1.write("average distribution of Ge-Pb (j around i) at 3rd shell: "+str(ave_pSRO_dist_ij_NN3)+'\n')
		f1.write("average distribution of Pb-Ge (i around j) at 3rd shell: "+str(ave_pSRO_dist_ji_NN3)+'\n')
	f1.close()
