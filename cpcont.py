import os
from shutil import copyfile

c=input("Sn concentration:")
data_folder = '/home/share/GeSn/64-atoms/2x2x2/Sn-{0}/MC-sampling'.format(c)

if not os.path.exists("CONTCAR_file"):
	os.mkdir("CONTCAR_file")

infile = open("process_nodoublecount.txt","r")
data = infile.readlines()
for line in data:
	index = str(line.split()[0])
	copyfile(os.path.join(data_folder,'{0}/CONTCAR'.format(index)),'./CONTCAR_file/CONTCAR.{0}'.format(index))
