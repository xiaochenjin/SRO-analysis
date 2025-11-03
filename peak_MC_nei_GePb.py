import os
import numpy as np
import matplotlib.pyplot as plt


def frange(start, end=None, inc=None):
#	"A range function, that does accept float increments..."

	if end == None:
		end = start + 0.0
		start = 0.0

	if inc == None:
		inc = 1.0

	L = []
	while 1:
		next = start + len(L) * inc
		if inc > 0 and next >= end:
			break
		elif inc < 0 and next <= end:
			break
		L.append(next)

	return L

conc_Pb_list = ['0.0625','0.09375','0.125']
N_exclude = 300
anal_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/analysis/'
data_dir = '/home/xcjin/Research/GePb/MC-sampling/128-atoms/data/'

if not os.path.exists(anal_dir):
	os.mkdir(anal_dir)

for conc_Pb in conc_Pb_list:
	print(conc_Pb)
	#N_div = int(input("Number of division: "))
	N_div=500

	#with open('./process_nodoublecount.txt') as f1:
	#with open('./energy_config.txt') as f1:
	#	index,energy = np.loadtxt(f1,dtype='str',delimiter =' ',usecols=(0,1),unpack = True)
	#index=['0000','0001','0002','0003','0004','0005']
	#index=['0000']

	data_conc = os.path.join(data_dir,conc_Pb)
	anal_conc = os.path.join(anal_dir,conc_Pb)
	if not os.path.exists(anal_conc):
		os.mkdir(anal_conc)

	os.chdir(data_conc)

	infile = open('./process_doublecount.txt','r')
	data=infile.readlines()[N_exclude:]
	index=[]
	for line in data:
		if "Accepted" in line:
			index.append(str(line.split()[0]))

	#lc = 5.745
	num_cell = 64

	neigh_dis=[np.sqrt(3/16),np.sqrt(1/2),np.sqrt(11/16),1,np.sqrt(19/16)]

	R_totallist=[];M_totallist=[];GR_totallist=[]
	num_bond_list=[];num_1_list=[];num_2_list=[];num_3_list=[];num_4_list=[];num_5_list=[]
	dist_per_num1=[];dist_per_num2=[];dist_per_num3=[];dist_per_num4=[];dist_per_num5=[]
	#dist_ave_list=[];dist_14_ave_list=[];dist_12_ave_list=[];dist_13_ave_list=[];dist_15_ave_list=[]
	#dist_1_ave_list=[];dist_2_ave_list=[];dist_3_ave_list=[];dist_4_ave_list=[];dist_5_ave_list=[]
	for cf in index:
		#READ NUMBER OF ATOMS FROM CONTCAR
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

		#Perfect location
		dist_per_num1.append(lc*neigh_dis[0])
		dist_per_num2.append(lc*neigh_dis[1])
		dist_per_num3.append(lc*neigh_dis[2])
		dist_per_num4.append(lc*neigh_dis[3])
		dist_per_num5.append(lc*neigh_dis[4])


		#Box length
		length = 2*lc
	#	print(length)
		dr = (length/2-0.01)/N_div
	#	dr = 0.06
		R_list=[]
		for R in frange(0.01,length/2,dr):
			R_list.append(R)
		R_totallist.append(R_list)


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

		
		#Calculation of RDF
		dist_list=[];M=N_div*[0];GR_list=[]
		density_0=N_Pb/V 
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
				k=int((rij-0.01)/dr) #k:index of R
			#	print(rij,k)
				if (k<N_div):
					M[k]=M[k]+1
		for l in range(N_div):
			density_R=(M[l]/N)/(4*3.14*R_list[l]*R_list[l]*dr)
			gr=density_R/density_0
			GR_list.append(gr)
		M_totallist.append(M)
		GR_totallist.append(GR_list)

	#	n_count=0
	#	for NR in range(1,N_div-1):
	#		if M[NR]<10 and M[NR]-M[NR-1]<0 and M[NR+1]-M[NR]>=0:
	#			n_count +=1
	#			if n_count ==1:
	#				NR_num1 = NR #upper boundary of 1st coordination shell)
	#				num_1 = np.sum(M[:NR_num1])
	#			if n_count ==2:
	#				NR_num2 = NR
	#				num_2 = np.sum(M[NR_num1:NR_num2])
	#			if n_count ==3:
	#				NR_num3 = NR
	#				num_3 = np.sum(M[NR_num2:NR_num3])

	#	num_bond_list.append(num_1/2)
	#	num_1_list.append(num_1)
	#	num_2_list.append(num_2)
	#	num_3_list.append(num_3)
	i
	#	print(num_1)
	#	print(num_2)
	#	print(num_3)


	#	plt.plot(R_list,M,'kx-')
	#	plt.title("Index: {0}".format(cf))
	#	plt.show()

	#CALCULATE AVERAGE
	R_list_ave = np.average(R_totallist,axis=0)
	M_ave = np.average(M_totallist,axis=0)
	GR_ave = np.average(GR_totallist,axis=0)
	#print(M_ave) 
	#print(len(M_ave))

	#Left boundary of peaks
	n_count=0
	for NR in range(1,N_div-1):
	#	if M_ave[NR]<15 and M_ave[NR]-M_ave[NR-1]<=0 and M_ave[NR+1]-M_ave[NR]>0:
		if M_ave[NR]-M_ave[NR-1]<=0 and M_ave[NR+1]-M_ave[NR]>0:
			n_count +=1
			if n_count ==1:
				NR_num1_left = NR #left boundary of 1st coordination shell)
			if n_count ==2:
				NR_num2_left = NR
			if n_count ==3:
				NR_num3_left = NR
			if n_count ==4:
				NR_num4_left = NR
			if n_count ==5:
				NR_num5_left = NR


	#Right boundary of peaks
	n_count=0
	for NR in range(1,N_div-1):
	#	if M_ave[NR]<15 and M_ave[NR]-M_ave[NR-1]<0 and M_ave[NR+1]-M_ave[NR]>=0:
		if M_ave[NR]-M_ave[NR-1]<0 and M_ave[NR+1]-M_ave[NR]>=0:
			n_count +=1
			if n_count ==1:
				NR_num1_right = NR #upper boundary of 1st coordination shell)
			if n_count ==2:
				NR_num2_right = NR
			if n_count ==3:
				NR_num3_right = NR
			if n_count ==4:
				NR_num4_right = NR
			if n_count ==5:
				NR_num5_right = NR


	#Calculate average NN and peak wide
	num_1_ave = np.sum(M_ave[NR_num1_left:NR_num1_right])
	num_2_ave = np.sum(M_ave[NR_num2_left:NR_num2_right])
	num_3_ave = np.sum(M_ave[NR_num3_left:NR_num3_right])

	try:
		NR_num4_left and NR_num4_right
		num_4_ave = np.sum(M_ave[NR_num4_left:NR_num4_right])
	except NameError:
		num_4_ave = 0

	try:
		NR_num5_left and NR_num5_right
		num_5_ave = np.sum(M_ave[NR_num5_left:NR_num5_right])
	except NameError:
		num_5_ave = 0


	width_peak_num1=round((R_list_ave[NR_num1_right]-R_list_ave[NR_num1_left])/2,3)
	width_peak_num2=round((R_list_ave[NR_num2_right]-R_list_ave[NR_num2_left])/2,3)
	width_peak_num3=round((R_list_ave[NR_num3_right]-R_list_ave[NR_num3_left])/2,3)

	try:
		NR_num4_left and NR_num4_right
		width_peak_num4=round((R_list_ave[NR_num4_right]-R_list_ave[NR_num4_left])/2,3)
	except NameError:
		width_peak_num4  = 0

	try:
		NR_num5_left and NR_num5_right
		width_peak_num5=round((R_list_ave[NR_num5_right]-R_list_ave[NR_num5_left])/2,3)
	except NameError:
		width_peak_num5 = 0


	dist_peak_num1=round((R_list_ave[NR_num1_right]+R_list_ave[NR_num1_left])/2,3)
	dist_peak_num2=round((R_list_ave[NR_num2_right]+R_list_ave[NR_num2_left])/2,3)
	dist_peak_num3=round((R_list_ave[NR_num3_right]+R_list_ave[NR_num3_left])/2,3)

	try:
		NR_num4_left and NR_num4_right
		dist_peak_num4=round((R_list_ave[NR_num4_right]+R_list_ave[NR_num4_left])/2,3)
	except NameError:
		dist_peak_num4  = 0

	try:
		NR_num5_left and NR_num5_right
		dist_peak_num5=round((R_list_ave[NR_num5_right]+R_list_ave[NR_num5_left])/2,3)
	except NameError:
		dist_peak_num5  = 0



	#Peak position
	##n_count=0
	#for NR in range(1,N_div-1):
	#	if M_ave[NR]-M_ave[NR-1]>0 and M_ave[NR+1]-M_ave[NR]<=0:
	#		n_count +=1
	#		if n_count ==1:
	#			dist_peak_num1=round(R_list_ave[NR],3)
	#		if n_count ==2:
	#			dist_peak_num2=round(R_list_ave[NR],3)
	#		if n_count ==3:
	#			dist_peak_num3=round(R_list_ave[NR],3)


	CN_GePb = round(num_1_ave/N_Pb,3)
	NN2_GePb = round(num_2_ave/N_Pb,3)
	NN3_GePb = round(num_3_ave/N_Pb,3)
	NN4_GePb = round(num_4_ave/N_Pb,3)
	NN5_GePb = round(num_5_ave/N_Pb,3)


	CN_PbGe = round(num_1_ave/N_Ge,3)
	NN2_PbGe = round(num_2_ave/N_Ge,3)
	NN3_PbGe = round(num_3_ave/N_Ge,3)
	NN4_PbGe = round(num_4_ave/N_Ge,3)
	NN5_PbGe = round(num_5_ave/N_Ge,3)


	#print(R_list_ave[NR_num1])
	#print(M_ave)

	perf_peak_num1 = round(np.average(dist_per_num1),3)
	perf_peak_num2 = round(np.average(dist_per_num2),3)
	perf_peak_num3 = round(np.average(dist_per_num3),3)
	perf_peak_num4 = round(np.average(dist_per_num4),3)
	perf_peak_num5 = round(np.average(dist_per_num5),3)

	#print("Ge around Pb Coordination number: ",CN_GePb)
	#print("Pb around Ge Coordination number: ",CN_PbGe)
	#print("perfect peak position CN: ", perf_peak_num1)
	#print("peak position CN: ",dist_peak_num1)
	#print("half peak width of 1st NN: ",width_peak_num1)
	#print("\n")
	#print("2nd Ge around Pb NN: ",NN2_GePb)
	#print("2nd Pb around Ge NN: ",NN2_PbGe)
	#print("perfect peak position NN2: ", perf_peak_num2)
	#print("peak position NN2: ",dist_peak_num2)
	#print("half peak width of 2nd NN: ",width_peak_num2)
	#print("\n")
	#print("3rd Ge around Pb NN: ",NN3_GePb)
	#print("3rd Pb around Ge NN: ",NN3_PbGe)
	#print("perfect peak position NN3: ", perf_peak_num3)
	#print("peak position NN3: ",dist_peak_num3)
	#print("half peak width of 3rd NN: ",width_peak_num3)


	num1_info = [CN_GePb,CN_PbGe,perf_peak_num1,dist_peak_num1,width_peak_num1]
	num2_info = [NN2_GePb,NN2_PbGe,perf_peak_num2,dist_peak_num2,width_peak_num2]
	num3_info = [NN3_GePb,NN3_PbGe,perf_peak_num3,dist_peak_num3,width_peak_num3]
	num4_info = [NN4_GePb,NN4_PbGe,perf_peak_num4,dist_peak_num4,width_peak_num4]
	num5_info = [NN5_GePb,NN5_PbGe,perf_peak_num5,dist_peak_num5,width_peak_num5]

	N_infor = len(num1_info)
	os.chdir(anal_conc)
	with open ("MC_peak_information_GePb.txt","w") as f1:
		for i in range(N_infor):
			f1.write(str(num1_info[i])+' '+str(num2_info[i])+' '+str(num3_info[i])+' '+str(num4_info[i])+' '+str(num5_info[i])+'\n')
	f1.close()

	N_seg = len(R_list_ave)
	with open ("MC_peak_NN_vs_distance_GePb.txt","w") as f1:
		for i in range(N_seg):
			f1.write(str(R_list_ave[i])+' '+str(M_ave[i])+'\n')
	f1.close()

	with open ("MC_RDF_NN_vs_distance_GePb.txt","w") as f1:
		for i in range(N_seg):
			f1.write(str(R_list_ave[i])+' '+str(GR_ave[i])+'\n')
	f1.close()



	#plt.xlabel('Distance (Angstrom) ')
	#plt.title('Total number of Ge-Pb NN vs Distance')
	#plt.ylabel('Total GePb NN')
	#plt.plot(R_list_ave,M_ave,'kx-')
	#plt.show()







	#with open("data_distance.txt","w") as f1:
	#	cf = 0
	#	for i in index:
	#		f1.write(str(i)+' '+str(dist_ave_list[cf])+' '+str(dist_12_ave_list[cf])+' '+str(dist_13_ave_list[cf])+' '+str(dist_14_ave_list[cf])+' '+str(dist_15_ave_list[cf])+'\n')
	#		cf=cf+1
	#f1.close()

	#with open("peak_data_neighbor_energy.txt","w") as f1:
	#	cf = 0
	#	for i in index:
	#		f1.write(str(i)+' '+str(num_bond_list[cf])+' '+str(num_1_list[cf])+' '+str(num_2_list[cf])+' '+str(num_3_list[cf])+' '+str(energy[cf])+'\n')
	#		cf=cf+1
	#f1.close()

	#with open("data_neighbor_distance.txt","w") as f1:
	#	cf = 0
	#	for i in index:
	#		f1.write(str(i)+' '+' '+str(dist_1_ave_list[cf])+' '+str(dist_2_ave_list[cf])+' '+str(dist_3_ave_list[cf])+' '+str(dist_4_ave_list[cf])+' '+str(dist_5_ave_list[cf])+'\n')
	#		cf=cf+1
	#f1.close()

	#dist_first = np.nanmean(dist_1_ave_list)
	#dist_second = np.nanmean(dist_2_ave_list)
	#dist_third = np.nanmean(dist_3_ave_list)
	#dist_forth = np.nanmean(dist_4_ave_list)
	#dist_fifth = np.nanmean(dist_5_ave_list)
	#
	#with open("neighbor_distance.txt","w") as f1:
	#	f1.write(str(dist_first)+' '+str(dist_second)+' '+str(dist_third)+' '+str(dist_forth)+' '+str(dist_fifth)+'\n')	
	#f1.close()
			
	#
