import numpy as np

#/Users/xcjin/WORK/GWU/Research/Paper/SRO in SiGeSn/Results/results-summary/average-basin/SRO-parameters
#conc_dict = {'Si': 0.25, "Ge":0.25,"Sn": 0.5}
#conc_dict = {'Si': 0.5, "Ge":0.5,"Sn": 1.0}
conc_dict = {'Si': 0.25, "Ge":0.5,"Sn": 0.25}

#ave_basin_pCN_*.txt
#second column (R-SRO) ; first column is E-SRO, last column is the average
KNN = '1'
if KNN == '1':
	pSRO_dict = {"SnSn": 0.66,"SiSi": -0.1,"GeGe": 0.17,"SiSn": 0.06,"SiGe": 0.02,"GeSn": -0.36} #SiSi and GeGe are the same
	CN_total = 4

CN_dict = {}

#around Sn
CN_dict['SiSn'] = CN_total*conc_dict['Si']*(1 - pSRO_dict['SiSn']) #Si around Sn
CN_dict['GeSn'] = CN_total*conc_dict['Ge']*(1 - pSRO_dict['GeSn']) 
CN_dict['SnSn'] = CN_total*conc_dict['Sn']*(1 - pSRO_dict['SnSn'])

print("Number of atoms around Sn: ", CN_dict['SiSn'] + CN_dict['GeSn'] + CN_dict['SnSn'])


#around Si
CN_dict['SiSi'] = CN_total*conc_dict['Si']*(1 - pSRO_dict['SiSn']) 
CN_dict['GeSi'] = CN_total*conc_dict['Ge']*(1 - pSRO_dict['SiGe'])
CN_dict['SnSi'] = CN_total*conc_dict['Sn']*(1 - pSRO_dict['SiSn'])

print("Number of atoms around Si: ", CN_dict['SiSi'] + CN_dict['GeSi'] + CN_dict['SnSi'])

#around Ge
CN_dict['SiGe'] = CN_total*conc_dict['Si']*(1 - pSRO_dict['SiGe'])
CN_dict['SnGe'] = CN_total*conc_dict['Sn']*(1 - pSRO_dict['GeSn'])
CN_dict['GeGe'] = CN_total*conc_dict['Ge']*(1 - pSRO_dict['GeGe'])

print("Number of atoms around Ge: ", CN_dict['SiGe'] + CN_dict['SnGe'] + CN_dict['GeGe'])


	
