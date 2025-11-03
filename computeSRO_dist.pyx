import os
from libc.math cimport sqrt
#from libc.math cimport fabs
cimport numpy as np
import numpy as np
#from cython.parallel import prange

neigh_dis=[sqrt(3/int(16)),sqrt(1/int(2)),sqrt(11/int(16)),1,sqrt(19/int(16))]

cpdef tuple PSRO(str n_f,int N_species,int N_shell,int num_cell,str type_pair,str data_dir,str conc_Pb, list conc_Pb_list,list BC_AA_totlist, list BC_AB_totlist): #list is the return type
#a cdef for C types and a def fr Python types (like numpy), cpdef is combined
	data_conc = os.path.join(data_dir,conc_Pb)
	infile = open(os.path.join(data_conc,'contcar-only/CONTCAR-{0}'.format(n_f)),'r')
	data = infile.readlines()
	infile.close()
	cdef list atomspecies = []
	cdef int N=0,i,j
	cdef np.ndarray H, H_T
	for i in range(N_species):
		atomspecies.append(int(data[6].split()[i]))
		N += int(data[6].split()[i])
	cdef float coe
	cdef list a=[],b=[],c=[]
	#READ COEFFICIENT
	coe = float(data[1].split()[0])
	for j in range(3):
#	for (int j=0;j<3;++j){
		a.append(float(data[2].split()[j])*coe)
		b.append(float(data[3].split()[j])*coe)
		c.append(float(data[4].split()[j])*coe)
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX

	#LATTICE CONSTANT
	cdef float V,lc
	V = np.dot(a,np.cross(b,c))
#	lc = (V/int(num_cell))**(1/int(3))
	lc = abs((4*V/int(num_cell))**(1/int(3)))
#	print(lc)


	#BASED ON xyz_POSCAR.py line 50
	cdef int n,m,N_atom,N_prev
	cdef list s,s_list,s_totlist=[]
	cdef np.ndarray s_vector
	for n in range(N_species):
		N_atom = atomspecies[n]
		s_list = [] #atomic positions for each species
		for i in range(N_atom):
			s = []
			for j in range(3):
				s.append(float(data[8+sum(atomspecies[:n])+i].split()[j]))
			s_vector = np.array([s[0],s[1],s[2]])
			s_list.append(s_vector) #atomic postions for one species
			#print(s_vector)
		s_totlist.append(s_list) #atomic positions for all species

	cdef int num_1,num_2,num_3
	cdef int ni,nj,Ni,Nj
	cdef float ci,cj #concentration of species i and j
	cdef float s_xij,s_yij,s_zij,rij
	cdef np.ndarray r_vector,si,sj
	cdef list si_list,sj_list
	cdef float CN1,CN2,CN3
	cdef float p1,p2,p3
	cdef list p_list,p_totlist=[],p_std_list,p_std_totlist=[]
	cdef list p1s_list1,p2s_list1,p3s_list1
	cdef list p1s_list2,p2s_list2,p3s_list2
	cdef list num1_ij_list,num2_ij_list,num3_ij_list
	cdef list num1_ji_list,num2_ji_list,num3_ji_list
	cdef list num1_ii_list,num2_ii_list,num3_ii_list
	cdef list num1_jj_list,num2_jj_list,num3_jj_list
	cdef list a1_ij,a1_ji,a2_ij,a2_ji,a3_ij,a3_ji #SRO parameter list
	cdef list a1_ii,a1_jj,a2_ii,a2_jj,a3_ii,a3_jj
	cdef list x1_ij,x1_ji,x2_ij,x2_ji,x3_ij,x3_ji #possible values of SRO parameters
	cdef list x1_ii,x1_jj,x2_ii,x2_jj,x3_ii,x3_jj
	cdef list y1_ij,y1_ji,y2_ij,y2_ji,y3_ij,y3_ji #corresponding counting
	cdef list y1_ii,y1_jj,y2_ii,y2_jj,y3_ii,y3_jj
	cdef list y_list,y_totlist = []
	cdef int index_pair

	if type_pair == 'A-B':
		index_pair = -1 #keep track of index of AB pairs
		for ni in range(N_species-1):
			Ni = atomspecies[ni]
			ci = Ni/int(N)
			si_list = s_totlist[ni]
			for nj in range(ni+1,N_species): 
				index_pair += 1
				BC_AB_CN = BC_AB_totlist[0][index_pair][conc_Pb_list.index(conc_Pb)]
				BC_AB_NN2 = BC_AB_totlist[1][index_pair][conc_Pb_list.index(conc_Pb)]
				BC_AB_NN3 = BC_AB_totlist[2][index_pair][conc_Pb_list.index(conc_Pb)]
				Nj = atomspecies[nj] 
				cj = Nj/int(N)
				sj_list = s_totlist[nj]
				num_1=0;num_2=0;num_3=0 #for one pair
				p1s_list1 = [];p2s_list1 = [];p3s_list1 = [] #list of SRO parameters for i atoms around each j atom
				num1_ji_list = [];num2_ji_list = [];num3_ji_list = [] #list of number of i atoms around each j atom
				for j in range(Nj):
					sj = sj_list[j]
					num_1s = 0;num_2s = 0;num_3s = 0 #number of i around each single j atom
					for i in range(Ni):
						si = si_list[i]
						s_xij=sj[0]-si[0]
						s_yij=sj[1]-si[1]
						s_zij=sj[2]-si[2]
						s_xij=s_xij-round(s_xij)
						s_yij=s_yij-round(s_yij)
						s_zij=s_zij-round(s_zij)
						s_vector=np.array([s_xij,s_yij,s_zij])
						r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
						rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])
#						if abs(rij-lc*neigh_dis[0])<0.5: num_1 +=1;num_1s += 1
#						if abs(rij-lc*neigh_dis[1])<0.35: num_2 +=1;num_2s += 1
#						if abs(rij-lc*neigh_dis[2])<0.35: num_3 +=1; num_3s += 1
						if BC_AB_CN[0]<= rij <=BC_AB_CN[1]:num_1 +=1;num_1s += 1
						if BC_AB_NN2[0]< rij <=BC_AB_NN2[1]: num_2 +=1;num_2s += 1
						if BC_AB_NN3[0]< rij <=BC_AB_NN3[1]: num_3 +=1; num_3s += 1

					num1_ji_list.append(num_1s);num2_ji_list.append(num_2s);num3_ji_list.append(num_3s)
					#SRO parameter for i atoms around each single j atom:
					p1s = 1 - num_1s/float(4*ci); p2s = 1 - num_2s/float(12*ci); p3s = 1 - num_3s/float(12*ci);
					p1s_list1.append(p1s);p2s_list1.append(p2s);p3s_list1.append(p3s)

				p1s_list2 = [];p2s_list2 = [];p3s_list2 = [] #list of SRO parameters for i atoms around each i atom
				num1_ij_list = [];num2_ij_list = [];num3_ij_list = [] #list of number of j atoms around each i atom
				for i in range(Ni):
					si = si_list[i]
					num_1s = 0;num_2s = 0;num_3s = 0 #number of j atoms around each single i atom
					for j in range(Nj):
						sj = sj_list[j]
						s_xij=sj[0]-si[0]
						s_yij=sj[1]-si[1]
						s_zij=sj[2]-si[2]
						s_xij=s_xij-round(s_xij)
						s_yij=s_yij-round(s_yij)
						s_zij=s_zij-round(s_zij)
						s_vector=np.array([s_xij,s_yij,s_zij])
						r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
						rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])
						#if abs(rij-lc*neigh_dis[0])<0.5: num_1s += 1
						#if abs(rij-lc*neigh_dis[1])<0.35: num_2s += 1
						#if abs(rij-lc*neigh_dis[2])<0.35: num_3s += 1
						if BC_AB_CN[0]<= rij <=BC_AB_CN[1]:num_1s += 1
						if BC_AB_NN2[0]< rij <=BC_AB_NN2[1]:num_2s += 1
						if BC_AB_NN3[0]< rij <=BC_AB_NN3[1]:num_3s += 1
					num1_ij_list.append(num_1s);num2_ij_list.append(num_2s);num3_ij_list.append(num_3s)
					#SRO parameter for j atoms around each single i atom:
					p1s = 1 - num_1s/float(4*cj); p2s = 1 - num_2s/float(12*cj); p3s = 1 - num_3s/float(12*cj);
					p1s_list2.append(p1s);p2s_list2.append(p2s);p3s_list2.append(p3s)

				CN1 = num_1/int(Nj); CN2 = num_2/int(Nj);CN3 = num_3/int(Nj) #average CNs for one pair ij
			#	p1 = CN1/4; p2 = CN2/12; p3 = CN3/12;
				p1 = 1 - CN1/float(4*ci); p2 = 1 - CN2/float(12*ci); p3 = 1 - CN3/float(12*ci);
				if N_shell == 1: p_list = [p1] #for one i-j pair
				if N_shell == 2: p_list = [p1,p2]
				if N_shell == 3: p_list = [p1,p2,p3]
				p_totlist.append(p_list) #for all i-j pairs			

				#standard deviation of SRO parameter for each single i atom:
				p1_std1 = np.std(p1s_list1);p2_std1 = np.std(p2s_list1);p3_std1 = np.std(p3s_list1);
				p1_std2 = np.std(p1s_list2);p2_std2 = np.std(p2s_list2);p3_std2 = np.std(p3s_list2);
				if N_shell == 1: p_std_list = [p1_std2,p1_std1] #for one i-j and j-i pair
				if N_shell == 2: p_std_list = [p1_std2,p1_std1,p2_std2,p2_std1]
				if N_shell == 3: p_std_list = [p1_std2,p1_std1,p2_std2,p2_std1,p3_std2,p3_std1]
				p_std_totlist.append(p_std_list) #for all i-j pairs	

				#compute distribution: a are SRO parameters 
				a1_ij = []; a1_ji=[];
				for i in range(len(num1_ij_list)): #list of number of j atoms around each i atom
					a1_ij.append(1 - num1_ij_list[i]/float(4*cj))
				for i in range(len(num1_ji_list)): #list of number of i atoms around eah j atom
					a1_ji.append(1 - num1_ji_list[i]/float(4*ci))

				a2_ij = [];a2_ji=[];
				for i in range(len(num2_ij_list)):
					a2_ij.append(1 - num2_ij_list[i]/float(12*cj))
				for i in range(len(num2_ji_list)):
					a2_ji.append(1 - num2_ji_list[i]/float(12*ci))

				a3_ij = []; a3_ji = []; 
				for i in range(len(num3_ij_list)):
					a3_ij.append(1 - num3_ij_list[i]/float(12*cj))
				for i in range(len(num3_ji_list)):
					a3_ji.append(1 - num3_ji_list[i]/float(12*ci))


				x1_ij = []
				for i in range(5):
					x1_ij.append(1 - i/float(4*cj))

				y1_ij = [0]*5
				for i in range(len(a1_ij)):
					for j in range(5):
						if a1_ij[i] == x1_ij[j]:
							y1_ij[j] += 1

				x2_ij = []
				for i in range(13):
					x2_ij.append(1 - i/float(12*cj))
			#   print(x2_ij)

				y2_ij = [0]*13
				for i in range(len(a2_ij)):
					for j in range(13):
						if a2_ij[i] == x2_ij[j]:
							y2_ij[j] += 1

				x3_ij = []
				for i in range(13):
					x3_ij.append(1 - i/float(12*cj))

				y3_ij = [0]*13
				for i in range(len(a3_ij)):
					for j in range(13):
						if a3_ij[i] == x3_ij[j]:
							y3_ij[j] += 1

				x1_ji = []
				for i in range(5):
						x1_ji.append(1 - i/float(4*ci))

				y1_ji = [0]*5
				for i in range(len(a1_ji)):
					for j in range(5):
						if a1_ji[i] == x1_ji[j]:
							y1_ji[j] += 1

				x2_ji = []
				for i in range(13):
					x2_ji.append(1 - i/float(12*ci))

				y2_ji = [0]*13
				for i in range(len(a2_ji)):
					for j in range(13):
						if a2_ji[i] == x2_ji[j]:
							y2_ji[j] += 1

				x3_ji = []
				for i in range(13):
					x3_ji.append(1 - i/float(12*ci))

				y3_ji = [0]*13
				for i in range(len(a3_ji)):
					for j in range(13):
						if a3_ji[i] == x3_ji[j]:
							y3_ji[j] += 1

				if N_shell == 1: y_list = [y1_ij,y1_ji] #for one i-j and j-i pair
				if N_shell == 2: y_list = [y1_ij,y1_ji,y2_ij,y2_ji]
				if N_shell == 3: y_list = [y1_ij,y1_ji,y2_ij,y2_ji,y3_ij,y3_ji]
				y_totlist.append(y_list) #for all i-j pairs

	if type_pair == 'A-A':
		for ni in range(N_species):
			BC_AA_CN = BC_AA_totlist[0][ni][conc_Pb_list.index(conc_Pb)]
			BC_AA_NN2 = BC_AA_totlist[1][ni][conc_Pb_list.index(conc_Pb)]
			BC_AA_NN3 = BC_AA_totlist[2][ni][conc_Pb_list.index(conc_Pb)]
			Ni = atomspecies[ni]
			ci = Ni/int(N)
			#print(Ni,N,N/int(Ni))
			si_list = s_totlist[ni]
			num_1=0;num_2=0;num_3=0 #for one pair
			p1s_list1 = [];p2s_list1 = [];p3s_list1 = [] #list of SRO parameters for i atoms around each i atom
			num1_ii_list = [];num2_ii_list = [];num3_ii_list = [] #list of number of i atoms around each i atom
			for i in range(Ni):
				si = si_list[i]
				num_1s = 0;num_2s = 0;num_3s = 0 #number of i around each single i atom
				for j in range(Ni):
					sj = si_list[j] #also from si_list because this is from one species
					s_xij=sj[0]-si[0]
					s_yij=sj[1]-si[1]
					s_zij=sj[2]-si[2]
					s_xij=s_xij-round(s_xij)
					s_yij=s_yij-round(s_yij)
					s_zij=s_zij-round(s_zij)
					s_vector=np.array([s_xij,s_yij,s_zij])
					r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
					rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])
#					print(abs(rij-(lc*neigh_dis[0])))
#					if abs(rij-lc*neigh_dis[0])<0.5: num_1 +=1;num_1s += 1
#					if abs(rij-lc*neigh_dis[1])<0.35: num_2 +=1;num_2s += 1
#					if abs(rij-lc*neigh_dis[2])<0.35: num_3 +=1;num_3s += 1
					if BC_AA_CN[0]<= rij <=BC_AA_CN[1]:num_1 +=1;num_1s += 1
					if BC_AA_NN2[0]< rij <=BC_AA_NN2[1]: num_2 +=1;num_2s += 1
					if BC_AA_NN3[0]< rij <=BC_AA_NN3[1]: num_3 +=1; num_3s += 1
				num1_ii_list.append(num_1s);num2_ii_list.append(num_2s);num3_ii_list.append(num_3s)
				#SRO parameter for i atoms around each single i atom:
				p1s = 1 - num_1s/float(4*ci); p2s = 1 - num_2s/float(12*ci); p3s = 1 - num_3s/float(12*ci);
				p1s_list1.append(p1s);p2s_list1.append(p2s);p3s_list1.append(p3s)


			CN1 = num_1/int(Ni); CN2 = num_2/int(Ni);CN3 = num_3/int(Ni)
#			print(2*num_1,2*num_2,2*num_3)
#			print(CN1,CN2,CN3)
			p1 = 1 - CN1/float(4*ci); p2 = 1 - CN2/float(12*ci); p3 = 1 - CN3/float(12*ci);
			if N_shell == 1: p_list = [p1] #for one i-i pair
			if N_shell == 2: p_list = [p1,p2]
			if N_shell == 3: p_list = [p1,p2,p3]
			p_totlist.append(p_list) #for all i-i pairs

			#standard deviation of SRO parameter for each single i atom:
			p1_std1 = np.std(p1s_list1);p2_std1 = np.std(p2s_list1);p3_std1 = np.std(p3s_list1);
			if N_shell == 1: p_std_list = [p1_std1] #for one i-i pair
			if N_shell == 2: p_std_list = [p1_std1,p2_std1]
			if N_shell == 3: p_std_list = [p1_std1,p2_std1,p3_std1]
			p_std_totlist.append(p_std_list) #for all i-j pairs

			#compute distribution
			a1_ii = []
			for i in range(len(num1_ii_list)): #list of number of i atoms around each i atom
				a1_ii.append(1 - num1_ii_list[i]/float(4*ci))

			a2_ii = []
			for i in range(len(num2_ii_list)):
				a2_ii.append(1 - num2_ii_list[i]/float(12*ci))

			a3_ii = []
			for i in range(len(num3_ii_list)):
				a3_ii.append(1 - num3_ii_list[i]/float(12*ci))


			x1_ii = []
			for i in range(5):
				x1_ii.append(1 - i/float(4*ci))

			y1_ii = [0]*5
			for i in range(len(a1_ii)):
				#print(a1_ii[i])
				for j in range(5):
					if a1_ii[i] == x1_ii[j]:
						y1_ii[j] += 1

			x2_ii = []
			for i in range(13):
				x2_ii.append(1 - i/float(12*ci))
	#   print(x2_ii)

			y2_ii = [0]*13
			for i in range(len(a2_ii)):
	#			print(a2_ii[i])
				for j in range(13):
					if a2_ii[i] == x2_ii[j]:
						y2_ii[j] += 1

			x3_ii = []
			for i in range(13):
				x3_ii.append(1 - i/float(12*ci))

			y3_ii = [0]*13
			for i in range(len(a3_ii)):
	#			 print(a3_ii[i])
				for j in range(13):
					if a3_ii[i] == x3_ii[j]:
						y3_ii[j] += 1

			if N_shell == 1: y_list = [y1_ii] #for one i-i pair
			if N_shell == 2: y_list = [y1_ii,y2_ii]
			if N_shell == 3: y_list = [y1_ii,y2_ii,y3_ii]
			y_totlist.append(y_list) #for all i-i pairs


#	return p_totlist
#	return list(np.array(p_totlist).flatten()),list(np.array(p_std_totlist).flatten()) #make 2D array into 1D
	return list(np.array(p_totlist).flatten()),list(np.array(p_std_totlist).flatten()),y_totlist#list(np.array(y_totlist).flatten())		

