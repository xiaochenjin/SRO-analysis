import numpy as np

#info of the target: 0.09375/SRO-1881
N = 128
conc_j = 0.09375 
conc_i = 1 - conc_j
#average/distriutions of SRO parameters computed by GePb-pSRO-AB/AA.py
pSRO_ave_ij_1NN  = -0.1034482792019844 #average i-j SRO 1NN
pSRO_ave_ij_2NN  = -0.07279693486590035 #average i-j SRO 2NN
pSRO_ave_ij_3NN  = 0.04980842911877398 #average i-j SRO 3NN
pSRO_ave_ii_1NN  = 0.010701540857553482
pSRO_ave_ii_2NN  = 0.009116131589377692 #-0.000396377727156505
pSRO_ave_ii_3NN  = 0.002774468855932355
pSRO_ave_jj_1NN  = 1.0
pSRO_ave_jj_2NN  = 0.7037037014961243
pSRO_ave_jj_3NN = -0.4814814329147339
SRO_ij_1NN = '70x44x2x0x0' #distribution of j around i at 1NN shell
SRO_ji_1NN = '0x0x0x0x12' #distribution of i around j at 1NN shell
SRO_ij_1NN = [float(x) for x in SRO_ij_1NN.split('x')]
SRO_ji_1NN = [float(x) for x in SRO_ji_1NN.split('x')]
SRO_ij_2NN = '29x47x27x13x0x0x0x0x0x0x0x0x0' 
SRO_ji_2NN = '0x0x0x0x0x0x0x0x0x0x0x4x8'
SRO_ij_2NN = [float(x) for x in SRO_ij_2NN.split('x')]
SRO_ji_2NN = [float(x) for x in SRO_ji_2NN.split('x')]
SRO_ij_3NN = '31x52x28x4x1x0x0x0x0x0x0x0x0'
SRO_ji_3NN = '0x0x0x0x0x0x0x0x0x1x7x3x1'
SRO_ij_3NN = [float(x) for x in SRO_ij_3NN.split('x')]
SRO_ji_3NN = [float(x) for x in SRO_ji_3NN.split('x')]
SRO_ii_1NN = '0x0x2x44x70'
SRO_ii_1NN = [float(x) for x in SRO_ii_1NN.split('x')]
#SRO_ii_2NN = '0x0x0x0x0x0x0x0x0x13x27x44x25'
SRO_ii_2NN = '0x0x0x0x0x0x0x0x0x15x25x47x29'
SRO_ii_2NN = [float(x) for x in SRO_ii_2NN.split('x')]
SRO_ii_3NN = '0x0x0x0x0x0x0x0x1x5x31x53x26'
SRO_ii_3NN = [float(x) for x in SRO_ii_3NN.split('x')]
SRO_jj_1NN = '12x0x0x0x0'
SRO_jj_1NN = [float(x) for x in SRO_jj_1NN.split('x')]
SRO_jj_2NN = '8x4x0x0x0x0x0x0x0x0x0x0x0'
SRO_jj_2NN = [float(x) for x in SRO_jj_2NN.split('x')]
SRO_jj_3NN = '1x3x7x1x0x0x0x0x0x0x0x0x0'
SRO_jj_3NN = [float(x) for x in SRO_jj_3NN.split('x')]

#print(SRO_ij_1NN)

#check if the distribution satisfy number of atoms.
Ni = N*conc_i
Nj = N*conc_j
Ni_dist_ij_1NN = np.sum(SRO_ij_1NN) #compute number of center atoms based on distributions
Nj_dist_ji_1NN = np.sum(SRO_ji_1NN)
Ni_dist_ii_1NN = np.sum(SRO_ii_1NN)
Ni_dist_jj_1NN = np.sum(SRO_jj_1NN)
Ni_dist_ij_2NN = np.sum(SRO_ij_2NN)
Nj_dist_ji_2NN = np.sum(SRO_ji_2NN)
Ni_dist_ii_2NN = np.sum(SRO_ii_2NN)
Ni_dist_jj_2NN = np.sum(SRO_jj_2NN)
Ni_dist_ij_3NN = np.sum(SRO_ij_3NN)
Nj_dist_ji_3NN = np.sum(SRO_ji_3NN)
Ni_dist_ii_3NN = np.sum(SRO_ii_3NN)
Ni_dist_jj_3NN = np.sum(SRO_jj_3NN)
print("check if the distribution satisfy number of atoms.")
print("actual number of i atoms: ", Ni)
print("number of i atoms computed by distribution of 1NN: ", Ni_dist_ij_1NN, Ni_dist_ii_1NN)
print("number of i atoms computed by distribution of 2NN: ", Ni_dist_ij_2NN, Ni_dist_ii_2NN)
print("number of i atoms computed by distribution of 3NN: ", Ni_dist_ij_3NN, Ni_dist_ii_3NN)
print("actual number of j atoms: ", Nj)
print("number of j atoms computed by distribution of 1NN: ", Nj_dist_ji_1NN, Ni_dist_jj_1NN)
print("number of j atoms computed by distribution of 2NN: ", Nj_dist_ji_2NN, Ni_dist_jj_2NN)
print("number of j atoms computed by distribution of 3NN: ", Nj_dist_ji_3NN, Ni_dist_jj_3NN)
print("\n")

#check if the distribution gives the consistent average as computed by cython code
NN1_list = np.arange(5) #possible occupation of 1st shell
pSRO_ij_1NN_list = [] #possible pSRO of 1st shell (j around i)
pSRO_ji_1NN_list = [] #possible pSRO of 1st shell (i around j)
pSRO_ii_1NN_list = [] 
pSRO_jj_1NN_list = []
NN2_list = np.arange(13)
pSRO_ij_2NN_list = [] 
pSRO_ji_2NN_list = [] 
pSRO_ii_2NN_list = []
pSRO_jj_2NN_list = []
pSRO_ij_3NN_list = []
NN3_list = np.arange(13)
pSRO_ji_3NN_list = []
pSRO_ii_3NN_list = []
pSRO_jj_3NN_list = []

for NN1 in NN1_list:
	pSRO_ij_1NN = 1 - NN1/(4*conc_j)
	pSRO_ij_1NN_list.append(pSRO_ij_1NN)
for NN1 in NN1_list:
	pSRO_ji_1NN = 1 - NN1/(4*conc_i)
	pSRO_ji_1NN_list.append(pSRO_ji_1NN)
for NN1 in NN1_list:
	pSRO_ii_1NN = 1 - NN1/(4*conc_i)
	pSRO_ii_1NN_list.append(pSRO_ii_1NN)
for NN1 in NN1_list:
	pSRO_jj_1NN = 1 - NN1/(4*conc_j)
	pSRO_jj_1NN_list.append(pSRO_jj_1NN)
#print(pSRO_ij_1NN_list)

for NN2 in NN2_list:
	pSRO_ij_2NN = 1 - NN2/(12*conc_j)
	pSRO_ij_2NN_list.append(pSRO_ij_2NN)
for NN2 in NN2_list:
	pSRO_ji_2NN = 1 - NN2/(12*conc_i)
	pSRO_ji_2NN_list.append(pSRO_ji_2NN)
for NN2 in NN2_list:
	pSRO_ii_2NN = 1 - NN2/(12*conc_i)
	pSRO_ii_2NN_list.append(pSRO_ii_2NN)
for NN2 in NN2_list:
	pSRO_jj_2NN = 1 - NN2/(12*conc_j)
	pSRO_jj_2NN_list.append(pSRO_jj_2NN)

for NN3 in NN3_list:
	pSRO_ij_3NN = 1 - NN3/(12*conc_j)
	pSRO_ij_3NN_list.append(pSRO_ij_3NN)
for NN3 in NN3_list:
	pSRO_ji_3NN = 1 - NN3/(12*conc_i)
	pSRO_ji_3NN_list.append(pSRO_ji_3NN)
for NN3 in NN3_list:
	pSRO_ii_3NN = 1 - NN3/(12*conc_i)
	pSRO_ii_3NN_list.append(pSRO_ii_3NN)
for NN3 in NN3_list:
	pSRO_jj_3NN = 1 - NN3/(12*conc_j)
	pSRO_jj_3NN_list.append(pSRO_jj_3NN)

#average SRO parameter computed by distributions
pSRO_ave_ij_1NN_dist = np.sum(np.multiply(pSRO_ij_1NN_list,SRO_ij_1NN))/np.sum(SRO_ij_1NN)
pSRO_ave_ji_1NN_dist = np.sum(np.multiply(pSRO_ji_1NN_list,SRO_ji_1NN))/np.sum(SRO_ji_1NN)
pSRO_ave_ii_1NN_dist = np.sum(np.multiply(pSRO_ii_1NN_list,SRO_ii_1NN))/np.sum(SRO_ii_1NN)
pSRO_ave_jj_1NN_dist = np.sum(np.multiply(pSRO_jj_1NN_list,SRO_jj_1NN))/np.sum(SRO_jj_1NN)

pSRO_ave_ij_2NN_dist = np.sum(np.multiply(pSRO_ij_2NN_list,SRO_ij_2NN))/np.sum(SRO_ij_2NN)
pSRO_ave_ji_2NN_dist = np.sum(np.multiply(pSRO_ji_2NN_list,SRO_ji_2NN))/np.sum(SRO_ji_2NN)
pSRO_ave_ii_2NN_dist = np.sum(np.multiply(pSRO_ii_2NN_list,SRO_ii_2NN))/np.sum(SRO_ii_2NN)
pSRO_ave_jj_2NN_dist = np.sum(np.multiply(pSRO_jj_2NN_list,SRO_jj_2NN))/np.sum(SRO_jj_2NN)

pSRO_ave_ij_3NN_dist = np.sum(np.multiply(pSRO_ij_3NN_list,SRO_ij_3NN))/np.sum(SRO_ij_3NN)
pSRO_ave_ji_3NN_dist = np.sum(np.multiply(pSRO_ji_3NN_list,SRO_ji_3NN))/np.sum(SRO_ji_3NN)
pSRO_ave_ii_3NN_dist = np.sum(np.multiply(pSRO_ii_3NN_list,SRO_ii_3NN))/np.sum(SRO_ii_3NN)
pSRO_ave_jj_3NN_dist = np.sum(np.multiply(pSRO_jj_3NN_list,SRO_jj_3NN))/np.sum(SRO_jj_3NN)

print("check if the distribution gives the consistent average as computed by cython code")
print("average i-j 1NN parameter computed by cython: ", pSRO_ave_ij_1NN)
print("average i-j 1NN parameter computed by j around i distribution", pSRO_ave_ij_1NN_dist)
print("average i-j 1NN parameter computed by i around j distribution", pSRO_ave_ji_1NN_dist) 
print("average i-i 1NN parameter computed by cython: ", pSRO_ave_ii_1NN)
print("average i-i 1NN parameter computed by i around i distribution", pSRO_ave_ii_1NN_dist)
print("average j-j 1NN parameter computed by cython: ", pSRO_ave_jj_1NN)
print("average j-j 1NN parameter computed by j around j distribution", pSRO_ave_jj_1NN_dist)

print("average i-j 2NN parameter computed by cython: ", pSRO_ave_ij_2NN)
print("average i-j 2NN parameter computed by j around i distribution", pSRO_ave_ij_2NN_dist)
print("average i-j 2NN parameter computed by i around j distribution", pSRO_ave_ji_2NN_dist)
print("average i-i 2NN parameter computed by cython: ", pSRO_ave_ii_2NN)
print("average i-i 2NN parameter computed by i around i distribution", pSRO_ave_ii_2NN_dist)
print("average j-j 2NN parameter computed by cython: ", pSRO_ave_jj_2NN)
print("average j-j 2NN parameter computed by j around j distribution", pSRO_ave_jj_2NN_dist)

print("average i-j 3NN parameter computed by cython: ", pSRO_ave_ij_3NN)
print("average i-j 3NN parameter computed by j around i distribution", pSRO_ave_ij_3NN_dist)
print("average i-j 3NN parameter computed by i around j distribution", pSRO_ave_ji_3NN_dist)
print("average i-i 3NN parameter computed by cython: ", pSRO_ave_ii_3NN)
print("average i-i 3NN parameter computed by i around i distribution", pSRO_ave_ii_3NN_dist)
print("average j-j 3NN parameter computed by cython: ", pSRO_ave_jj_3NN)
print("average j-j 3NN parameter computed by j around j distribution", pSRO_ave_jj_3NN_dist)
print("\n")

print("Check if A-B/B-A distribution can give A-A/B-B distribution")
#If it holds, then sum of CN around A or B should hold
print("Distribution of i-i (i around i) 1NN SRO parameter computed directly: ", SRO_ii_1NN)
print("Distribution of i-i (i around i) 1NN SRO parameter computed from i-j (j around i) distribution: ", np.flip(SRO_ij_1NN)) 
#np.flip doesn't change orginal array
print("Distribution of j-j (j around j) 1NN SRO parameter computed directly: ", SRO_jj_1NN)
print("Distribution of j-j (j around j) 1NN SRO parameter computed from j-i (i around j) distribution: ", np.flip(SRO_ji_1NN))

print("Distribution of i-i (i around i) 2NN SRO parameter computed directly: ", SRO_ii_2NN)
print("Distribution of i-i (i around i) 2NN SRO parameter computed from i-j (j around i) distribution: ", np.flip(SRO_ij_2NN))
print("Distribution of j-j (j around j) 2NN SRO parameter computed directly: ", SRO_jj_2NN)
print("Distribution of j-j (j around j) 2NN SRO parameter computed from j-i (i around j) distribution: ", np.flip(SRO_ji_2NN))

print("Distribution of i-i (i around i) 3NN SRO parameter computed directly: ", SRO_ii_3NN)
print("Distribution of i-i (i around i) 3NN SRO parameter computed from i-j (j around i) distribution: ", np.flip(SRO_ij_3NN))
print("Distribution of j-j (j around j) 3NN SRO parameter computed directly: ", SRO_jj_3NN)
print("Distribution of j-j (j around j) 3NN SRO parameter computed from j-i (i around j) distribution: ", np.flip(SRO_ji_3NN))
print("\n")

print("Check if the average consistent with CN computed by peak-distance-*.py previously")
#data_neighbor_energy_PbPb.txt
#1881 0 0 4 20 8 12 -641.31852329
total_num_jj = 0 #total number of j-j 1NN (total number of j-j bond = half of this)
CN1_jj = total_num_jj/Nj #average number of j around j
CN1_ji = 4 - CN1_jj #average number of i around j
pSRO_ji_original = 1 - CN1_ji/(4*conc_i) #pSRO based on average computed originally
print("average i-j 1NN parameter computed based on original calculation of CN: ", pSRO_ji_original)
print("average i-j 1NN parameter computed by cython: ", pSRO_ji_1NN)
print("\n")



