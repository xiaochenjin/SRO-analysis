CF = 1460
start_line = 18

infile = open("log","r")
data = infile.readlines()

index_list=[];energy_list=[];
P_list=[]; result_list=[];


for cf in range(CF):
	line_index = data[start_line+4*cf]
	line_energy = data[start_line+1+4*cf]
#	print(line_index)
#	print(line_energy)
	line_proba = data[start_line+2+4*cf]
	line_result = data[start_line+3+4*cf]
	
	index = str(line_index.split()[1])[1:]
	index_list.append(index)
	energy = str(line_energy.split()[2])
	energy_list.append(energy)
	proba = str(line_proba.split()[2])
	P_list.append(proba)
	result = str(line_result.split()[0])
	result_list.append(result)
infile.close()

with open ("process_nodoublecount.txt","w") as f1:
	for cf in range(CF):
		f1.write(str(index_list[cf])+' '+str(energy_list[cf])+' '+str(P_list[cf])+' '+str(result_list[cf])+'\n')
f1.close()

infor_list = []
for cf in range(CF):
	infor = [index_list[cf],energy_list[cf],P_list[cf],result_list[cf]]
	if infor[3] == "Accepted":
		infor_list.append(infor)
	if infor[3] == "Declined":
		infor_list.append(infor)
		infor_list.append(infor_list[-2])

total_step = len(infor_list)
with open ("process_doublecount.txt","w") as f1:
	for step in range(total_step):
		infor = infor_list[step]
		step_index = int(infor[0])
		f1.write(str(infor[0])+' '+str(infor[1])+' '+str(infor[2])+' '+str(infor[3])+'\n')
f1.close()

