import numpy as np
import random
import pure_rdf
import read_xyz
from optparse import OptionParser

#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--dimension', type = int,default = 10,help = 'unit is nm (default: %default)')

(options, args) = parser.parse_args()
dimension = options.dimension #nm

#basic structure info
#species = ['Ge']

input_name = '{0}x{0}x{0}-cube-GeSn-raw-APT.xyz'.format(str(dimension))
#input_name = 'denser-cube-Ge-raw-APT.xyz'
frac_position_list,cell_geometry = read_xyz.read(input_name)
#print(cell_geometry)


N_bin = int(1000*dimension)
R_list,GR_list = pure_rdf.rdf(frac_position_list,cell_geometry,N_bin)

output_name = 'rdf-{0}x{0}x{0}-cube-GeSn.txt'.format(str(dimension))
#output_name = 'rdf-denser-cube-Ge.txt'
with open (output_name,"w") as f1:
	for i in range(N_bin):
		f1.write(str(R_list[i])+' '+str(GR_list[i])+'\n')
f1.close()




