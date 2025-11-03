import operator

#count occurance of each configurations:
infile = open ("process_doublecount.txt","r")
data = infile.readlines()

index_list=[]
energy_list=[]
for line in data:
	if "Accepted" in line:
		index =line.split()[0]
		index_list.append(index)
		energy = line.split()[1]
		energy_list.append(energy)
infile.close()

index_energy_list = list(zip(index_list,energy_list))
#print(index_energy_list)

#https://stackoverflow.com/questions/23240969/python-count-repeated-elements-in-the-list/23240989

#count_dic = {index:index_list.count(index) for index in index_list}
count_dic = {index_energy:index_energy_list.count(index_energy) for index_energy in index_energy_list}

#https://www.tutorialspoint.com/How-to-convert-Python-Dictionary-to-a-list
count_list=count_dic.items()
#print(count_list)

#https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
#https://www.geeksforgeeks.org/python-program-to-sort-a-list-of-tuples-by-second-item/

sorted_count_list = sorted(count_list, key=operator.itemgetter(1),reverse = True)
num_accept = len(count_list) #number of configurations have been at least accepted once

#print(sorted_count_list)

with open ("count_accept_config.txt","w") as f1:
	for i in range(num_accept):
		f1.write(str(sorted_count_list[i][0][0])+' '+str(sorted_count_list[i][0][1])+' '+str(sorted_count_list[i][1])+'\n')
f1.close()

infile = open("process_nodoublecount.txt","r")
data = infile.readlines()
index_reject=[]
energy_reject=[]
for line in data:
	if "Declined" in line:
		index_reject.append(str(line.split()[0]))
		energy_reject.append(str(line.split()[1]))
num_reject=len(index_reject)
infile.close()

with open ("count_accept_config.txt","a") as f1:
	for i in range(num_reject):
		f1.write(str(index_reject[i])+' '+str(energy_reject[i])+' '+'0'+'\n')
f1.close()

infile = open("count_accept_config.txt","r")
data = infile.readlines()
index_list=[]
count_list=[]
energy_list=[]
for line in data:
	index_list.append(str(line.split()[0]))
	energy_list.append(str(line.split()[1]))
	count_list.append(str(line.split()[2]))
infile.close()
	
index_energy_count = list(zip(index_list,energy_list,count_list))
sorted_index_energy_count = sorted(index_energy_count, key=operator.itemgetter(1),reverse = True)

num_config = len(sorted_index_energy_count)

with open ("energy_config.txt","w") as f1:
    for i in range(num_config):
        f1.write(str(sorted_index_energy_count[i][0])+' '+str(sorted_index_energy_count[i][1])+' '+str(sorted_index_energy_count[i][2])+'\n')
f1.close()




