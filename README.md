# trans.py
Convert MC trajectory log file to process_nodoublecount.txt and process_doublecount.txt 

# acceptcount.py
Count number of acceptance for each trajectory and sort the configuration based on the total energy (output: energy_config.txt)

# peak_MC_nei_*.py:
Compute peaks (number of bond counts,rdf..) vs distance 

# GePb-SRO-AA.py / GePb-SRO-AB.py
Compute SRO parameter for each accepted snapshot in a Monte Carlo trajectory (cindluding double counting of the accepted configuration)

# average-pSRO-AA.py / average-pSRO-AB.py
Compute ensemble-averaged SRO parameter

# check-SRO.py
Check if the SRO are computed correctly through checking wether A-A/B-B SRO parameters consistent with A-B SRO parameters

# check-error.py
Check if the average SRO parameter and SRO distribution are computed consistently.

# compute-all-rdf.py and pure_rdf.pyx
Compute all-atom RDF

# obtain_kNN_location.py 
Obtain KNN peak positions from RDF (output: kNN-peak-position.txt)

# peak_distance_*.py
Compute coordination numbers (up to 5NN) for each configuratino

# cpcont.py
Copy data

# plot-ave-SRO.py
Plot average SRO parameters vs composition

# plot-ave-AA-every-step.py
Plot average SRO parameter vs MC step

# plot-random-MC-dist.py
Plot SRO parameter distribution (random vs SRO)

# plot-theory-dist-*.py
Plot theorectical SRO distribution (random alloy)

# plot_*NN_dist.py
Plot KNN SRO distribution for SRO alloys modeled by MC sampling

# random-select.py
Random select configurations from MC trajectories

# setup.py
Compile *.pyx into *c and *.so
python setup.py build_ext --inplace



