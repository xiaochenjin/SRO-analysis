c_Pb = 0.09375
c_Ge = 1 - c_Pb

pSRO_1NN_GePb = -0.098
pSRO_2NN_GePb = -0.066
pSRO_3NN_GePb = 0.055

#Ge around Pb
CN1_PbGe = 4*c_Ge*(1-pSRO_1NN_GePb) 
CN2_PbGe = 12*c_Ge*(1-pSRO_2NN_GePb)
CN3_PbGe = 12*c_Ge*(1-pSRO_3NN_GePb)

CN1_PbPb = 4 - CN1_PbGe
CN2_PbPb = 12 - CN2_PbGe
CN3_PbPb = 12 - CN3_PbGe

pSRO_1NN_PbPb = 1 - CN1_PbPb/(4*c_Pb)
pSRO_2NN_PbPb = 1 - CN2_PbPb/(12*c_Pb)
pSRO_3NN_PbPb = 1 - CN3_PbPb/(12*c_Pb)


#Pb around Ge
CN1_GePb = 4*c_Pb*(1-pSRO_1NN_GePb)
CN2_GePb = 12*c_Pb*(1-pSRO_2NN_GePb)
CN3_GePb = 12*c_Pb*(1-pSRO_3NN_GePb)

CN1_GeGe = 4 - CN1_GePb
CN2_GeGe = 12 - CN2_GePb
CN3_GeGe = 12 - CN3_GePb

pSRO_1NN_GeGe = 1 - CN1_GeGe/(4*c_Ge)
pSRO_2NN_GeGe = 1 - CN2_GeGe/(12*c_Ge)
pSRO_3NN_GeGe = 1 - CN3_GeGe/(12*c_Ge)

print("average 1NN Ge-Ge SRO parameter: ", pSRO_1NN_GeGe)
print("average 2NN Ge-Ge SRO parameter: ", pSRO_2NN_GeGe)
print("average 3NN Ge-Ge SRO parameter: ", pSRO_3NN_GeGe)

print("average 1NN Pb-Pb SRO parameter: ", pSRO_1NN_PbPb)
print("average 2NN Pb-Pb SRO parameter: ", pSRO_2NN_PbPb)
print("average 3NN Pb-Pb SRO parameter: ", pSRO_3NN_PbPb)
