import os
import numpy as np

conc_Pb_list = ['0.0625','0.09375','0.125']
BC_CN_list=[[2.59,2.67],[2.59,2.69],[2.61,2.70]]
BC_NN2_list=[[3.88,4.20],[3.81,4.25],[3.83,4.31]]
BC_NN3_list=[[4.26,4.92],[4.26,5.06],[4.31,5.23]]

anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'
data_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/data/'

for conc_Pb in conc_Pb_list:
	print(conc_Pb)
	data_conc = os.path.join(data_dir,conc_Pb)
	anal_conc = os.path.join(anal_dir,conc_Pb)

	os.chdir(data_conc)

	#with open('./process_nodoublecount.txt') as f1:
	with open('./energy_config.txt') as f1:
		index,energy = np.loadtxt(f1,dtype='str',delimiter =' ',usecols=(0,1),unpack = True)

	#with open ('./peak_information_GePb.txt') as f2:
	#	num1_info,num2_info,num3_info=np.loadtxt(f2,delimiter =' ',usecols=(0,1,2),unpack = True)

	#peak_dist=num2_info[3] #peak position of 2nd NN
	#width=num2_info[4] # half of peak width of 2nd NN
	#peak_dist=(num2_info[3]+num3_info[3])/2 #peak position of 2nd NN
	#width=num2_info[4]+num3_info[4] # half of peak width of 2nd NN

	#print(peak_dist)
	#print(width)
	BC_CN = BC_CN_list[conc_Pb_list.index(conc_Pb)]
	BC_NN2 = BC_NN2_list[conc_Pb_list.index(conc_Pb)]
	BC_NN3 = BC_NN3_list[conc_Pb_list.index(conc_Pb)]


	#infile = open('./process_doublecount.txt','r')
	#data=infile.readlines()
	#index=[]
	#for line in data:
	#	if "Accepted" in line:
	#		index.append(str(line.split()[0]))

	#lc = 5.745
	num_cell = 64


	num_bond_list=[];num_1_list=[];num_2_list=[];num_3_list=[];num_4_list=[];num_5_list=[]
	dist_ave_list=[];dist_14_ave_list=[];dist_12_ave_list=[];dist_13_ave_list=[];dist_15_ave_list=[]
	dist_1_ave_list=[];dist_2_ave_list=[];dist_3_ave_list=[];dist_4_ave_list=[];dist_5_ave_list=[]
	for cf in index:
		#READ NUMBER OF ATOMS FROM CONTCAR
		if not os.path.exists('./contcar-only/CONTCAR-{0}'.format(cf)):
			energy.remove(energy[index.index(cf)]);index.remove(cf)
			continue
		infile= open('./contcar-only/CONTCAR-{0}'.format(cf), 'r') #CONTCAR
		data = infile.readlines()
		N_Ge = int(data[6].split()[0])
		N_Pb = int(data[6].split()[1])
		N = N_Ge + N_Pb
		#READ COEFFICIENT
		coe = float(data[1].split()[0])
		#REWRITE CELL MATRIX
		a = [] #3 CELL VECTORS
		b = []
		c = []
		for j in range(3):
			a.append(float(data[2].split()[j])*coe)
			b.append(float(data[3].split()[j])*coe)
			c.append(float(data[4].split()[j])*coe)
		H_T = np.array([a,b,c])
		H = H_T.transpose() #CELL MATRIX

		#LATTICE CONSTANT
		V = np.dot(a,np.cross(b,c))	
	#	lc = (V/num_cell)**(1/3)
		lc = abs((4*V/num_cell)**(1/3))
	#	print(lc)

		#BASED ON xyz_POSCAR.py line 50
		sx_1=[];sy_1=[];sz_1=[]
		for i in range(N_Ge):
			s = []
			for j in range(3):
				s.append(float(data[8+i].split()[j]))
			s_vector = np.array([s[0],s[1],s[2]])
			sx_1.append(s_vector[0])
			sy_1.append(s_vector[1])
			sz_1.append(s_vector[2])

		#BASED ON xyz_POSCAR.py line 50
		sx_2=[];sy_2=[];sz_2=[]
		for i in range(N_Pb):
			s = []
			for j in range(3):
				s.append(float(data[8+N_Ge+i].split()[j]))
			s_vector = np.array([s[0],s[1],s[2]])
			sx_2.append(s_vector[0])
			sy_2.append(s_vector[1])
			sz_2.append(s_vector[2])
		
		dist_list=[]
		for i in range(N_Ge):
			for j in range(N_Pb):
				s_xij=sx_2[j]-sx_1[i]
				s_yij=sy_2[j]-sy_1[i]
				s_zij=sz_2[j]-sz_1[i]
				s_xij=s_xij-round(s_xij)
				s_yij=s_yij-round(s_yij)
				s_zij=s_zij-round(s_zij)
				s_vector=np.array([s_xij,s_yij,s_zij])
				r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
				rij=np.sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])
				sij = np.sqrt(s_vector[0]*s_vector[0]+s_vector[1]*s_vector[1]+s_vector[2]*s_vector[2])
				dist_list.append(rij)
		dist_list2 = sorted(dist_list)
	#	print(dist_list2)
		N_distance = np.size(dist_list2)
	#	for l in range(N_distance-1):
	#		if(dist_list2[l+1]-dist_list2[l]>0.5):
	#			print(dist_list2[l]) #print upper limit of nst neighbors

		neigh_dis=[np.sqrt(3/16),np.sqrt(1/2),np.sqrt(11/16),1,np.sqrt(19/16)] #neighbor distance of perfect diamond cubic with unit lattice constant
		dist_1=[];dist_2=[];dist_3=[];dist_4=[];dist_5=[];
		for l in range(N_distance):	
	#	   if abs(dist_list[l]-lc*neigh_dis[0]-0.2)<0.3:
	#		if abs(dist_list[l]-lc*neigh_dis[0])<0.5:
			if BC_CN[0]<= dist_list[l]< BC_CN[1]:
				dist_1.append(dist_list[l])
	#		if abs(dist_list[l]-peak_dist)<width:
			if BC_NN2[0]<=dist_list[l]< BC_NN2[1]:
				dist_2.append(dist_list[l])
	#		if abs(dist_list[l]-lc*neigh_dis[2])<0.5:
			if BC_NN3[0]<=dist_list[l]< BC_NN3[1]:
				dist_3.append(dist_list[l])
			if abs(dist_list[l]-lc*neigh_dis[3])<0.5:
				dist_4.append(dist_list[l])
			if abs(dist_list[l]-lc*neigh_dis[4])<0.5:
				dist_5.append(dist_list[l])

		#number of nearest neighbors: double count
		num_1=np.size(dist_1);num_2=np.size(dist_2);num_3=np.size(dist_3);num_4=np.size(dist_4);num_5=np.size(dist_5)
		num_1_list.append(num_1);num_2_list.append(num_2);num_3_list.append(num_3);num_4_list.append(num_4);num_5_list.append(num_5)
		num_bond=int(np.size(dist_1)/2);num_bond_list.append(num_bond)
		dist_12 = dist_1+dist_2
		dist_13 = dist_1+dist_2+dist_3
		dist_14 = dist_1+dist_2+dist_3+dist_4
		dist_15 = dist_1+dist_2+dist_3+dist_4+dist_5
		dist_12_ave = np.average(dist_12);dist_12_ave_list.append(dist_12_ave)
		dist_13_ave = np.average(dist_13);dist_13_ave_list.append(dist_13_ave)
		dist_14_ave = np.average(dist_14);dist_14_ave_list.append(dist_14_ave)
		dist_15_ave = np.average(dist_15);dist_15_ave_list.append(dist_15_ave)
		dist_ave = np.average(dist_list);dist_ave_list.append(dist_ave)

		if dist_1: 
			dist_1_ave = np.average(dist_1);dist_1_ave_list.append(round(dist_1_ave,3))
		else:dist_1_ave_list.append(np.nan)

		dist_2_ave = np.average(dist_2);dist_2_ave_list.append(round(dist_2_ave,3))
		dist_3_ave = np.average(dist_3);dist_3_ave_list.append(round(dist_3_ave,3))
		dist_4_ave = np.average(dist_4);dist_4_ave_list.append(round(dist_4_ave,3))
		dist_5_ave = np.average(dist_5);dist_5_ave_list.append(round(dist_5_ave,3))
		
	#	print(dist_1)
	#	print(dist_2)
	#	print(dist_3)
	#	print(dist_4)
	#	print(dist_5)


	#with open("data_distance.txt","w") as f1:
	#	cf = 0
	#	for i in index:
	#		f1.write(str(i)+' '+str(dist_ave_list[cf])+' '+str(dist_12_ave_list[cf])+' '+str(dist_13_ave_list[cf])+' '+str(dist_14_ave_list[cf])+' '+str(dist_15_ave_list[cf])+'\n')
	#		cf=cf+1
	#f1.close()

	os.chdir(anal_conc)
	with open("data_neighbor_energy_GePb.txt","w") as f1:
		cf = 0
		for i in index:
			f1.write(str(i)+' '+str(num_bond_list[cf])+' '+str(num_1_list[cf])+' '+str(num_2_list[cf])+' '+str(num_3_list[cf])+' '+str(num_4_list[cf])+' '+str(num_5_list[cf])+' '+str(energy[cf])+' '+'\n')
			cf=cf+1
	f1.close()

	with open("data_neighbor_distance_GePb.txt","w") as f1:
		cf = 0
		for i in index:
			f1.write(str(i)+' '+' '+str(dist_1_ave_list[cf])+' '+str(dist_2_ave_list[cf])+' '+str(dist_3_ave_list[cf])+' '+str(dist_4_ave_list[cf])+' '+str(dist_5_ave_list[cf])+'\n')
			cf=cf+1
	f1.close()

	dist_first = np.nanmean(dist_1_ave_list)
	dist_second = np.nanmean(dist_2_ave_list)
	dist_third = np.nanmean(dist_3_ave_list)
	dist_forth = np.nanmean(dist_4_ave_list)
	dist_fifth = np.nanmean(dist_5_ave_list)

	with open("neighbor_distance_GePb.txt","w") as f1:
		f1.write(str(dist_first)+' '+str(dist_second)+' '+str(dist_third)+' '+str(dist_forth)+' '+str(dist_fifth)+'\n')	
	f1.close()
		

