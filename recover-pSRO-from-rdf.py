import os
import numpy as np

def recover(R_list,GR_list,n_shell):
#	N_total = 46656
#	a = 98.69202
#	V = a**3
#	if pair == 'SiSi': c_Sn = 0.7
#	if pair == 'GeGe': c_Sn = 0.3
    c_Sn = 1 #fictious concentration
	N_Sn = N_total*c_Sn
	density_0 = N_Sn/V
#	N_bin = 50000
#	dr = 0.5*a/N_bin
	dr = R_list[1] - R_list[0]
#	print(R_list[:10])

	CN_SnSn = 0
	for i in range(len(R_list)):
		gr = GR_list[i]
		R = R_list[i]
		density_R = density_0*gr
		CN_dr = density_R*(4*3.14*R*R*dr)
		CN_SnSn += CN_dr

	Nk = Nk_list[n_shell]
	pSRO_SnSn = 1 - CN_SnSn/(Nk*c_Sn) 

	return pSRO_SnSn

def compute_pSRO(R_totlist,GR_totlist):
	index_range_list = [] #find index of R's closes to R_range
	for n_shell in range(len(Nk_list)):
		if 0 < n_shell < len(Nk_list)-1:
			length_lower = 0.5*lc*(neigh_dis[n_shell]-neigh_dis[n_shell-1])
			length_upper = 0.5*lc*(neigh_dis[n_shell+1]-neigh_dis[n_shell])
			length = min(length_lower,length_upper)
		if n_shell == 0:
			length  = 0.5*lc*(neigh_dis[n_shell+1]-neigh_dis[n_shell])
		if n_shell == len(Nk_list)-1:
			length = 0.5*lc*(neigh_dis[n_shell]-neigh_dis[n_shell-1])
		R_upper = lc*neigh_dis[n_shell]+length
		R_lower = lc*neigh_dis[n_shell]-length
		R_range = [R_lower,R_upper]
		#https://stackoverflow.com/questions/9706041/finding-index-of-an-item-closest-to-the-value-in-a-list-thats-not-entirely-sort
		index_lower = min(range(len(R_totlist)), key=lambda i: abs(R_totlist[i]-R_lower))
		index_higher = min(range(len(R_totlist)), key=lambda i: abs(R_totlist[i]-R_upper))
		index_range = [index_lower,index_higher]
		index_range_list.append(index_range)

	pSRO_list = []
	for n_shell in range(len(index_range_list)):
		index_range = index_range_list[n_shell]
		R_list = R_totlist[index_range[0]:index_range[1]]
		GR_list = GR_totlist[index_range[0]:index_range[1]]
		pSRO = recover(R_list,GR_list,n_shell)
		pSRO_list.append(pSRO)

	return pSRO_list			

def grep_GR(folder):
	file_name = os.path.join(folder,"rdf-{0}-noperturb.txt".format(pair))
	with open (file_name) as f1:
            R_totlist,GR_totlist = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)
	return R_totlist,GR_totlist

with open ("/home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/optimization-test/obtain-kNN-infor/kNN-peak-position.txt") as f1:
	neigh_dis,Nk_list = np.loadtxt(f1, delimiter=' ', usecols=(1,2), unpack=True) #Nk is total number of nearest neighbor at shell K

pair = input("pair: ")
#folder = '../../'
folder = '../../analysis/recover-test-4/'

#from obtain-kNN-infor/obtain-lattice-constant.py
N_total = 46656
a = 98.69202
V = a**3
num_cell = 18
lc = a/num_cell

R_totlist_SRO,GR_totlist_SRO = grep_GR(folder)
pSRO_list_SRO = compute_pSRO(R_totlist_SRO,GR_totlist_SRO)

print(pSRO_list_SRO[:4])

with open (os.path.join(folder,"pSRO-{0}-noperturb.txt".format(pair)),'w') as f1:
	for i in range(len(pSRO_list_SRO)):
		f1.write(str(i+1)+'NN '+str(pSRO_list_SRO[i])+'\n')
f1.close()



	









