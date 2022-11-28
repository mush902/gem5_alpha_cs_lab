import os
#import progressbar
import time
import subprocess

import re

from os.path import exists as file_exists

import numpy as np

ITERATIONS = 1

return_code = 0

#synthetic = ['tornado', 'bit_complement', 'bit_reverse', 'bit_rotation', 'neighbor', 'shuffle', 'transpose']
synthetic = ['transpose']

#global_threshold = 0

for traffic in synthetic:
    if file_exists('Q_Table.txt'):
        print('Q_table exists')
        os_command = "rm Q_Table.txt"
        os.system(os_command)

    global_threshold = 0.0
    for injrate in np.arange(0.05, 1.0, 0.05):

        local_threshold = 100000.0
        #if file_exists('Q_Table.txt'):
        #  print('Q_table exists')
       
        #os_command = "rm Q_Table.txt"
        #os.system(os_command)
    
        #if not file_exists('Q_Table.txt'):
        #   print('Q_table deleted')
           
        if return_code == 0:
            print("Starting Training Phase")
        else:
	    continue
	
        for i in range(ITERATIONS):
            print("Iteration Number:", i)
            os_command = "./build/ALPHA/gem5.opt -d m5out/ configs/example/garnet_synth_traffic.py --num-cpus=64 --num-dirs=64 --network=garnet2.0 --topology=Mesh_XY --mesh-rows=8 --sim-cycles=10000 --vcs-per-vnet=8 --injectionrate={:f} --synthetic={:s}  --routing-algorithm=2".format(injrate, traffic)
	    result = os.system(os_command)

	    if file_exists('file.txt'):

    	        os_command = 'rm file.txt'

    	        os.system(os_command)

	    os_command = 'grep "average_packet_latency" m5out/stats.txt > file.txt'

	    os.system(os_command)

	    fp = open('file.txt' , 'r')

	    lines = fp.readlines()

	    for line in lines:
	        elem = re.findall("\d+\.\d+", line)
	        for val in elem:
	    	    print("Current value is : ", val)
		    #print(val)
	            #if float(val) > float(global_threshold) and float(val)<=float(local_threshold):
                    if float(val)<=float(local_threshold):
            	        local_threshold = val
		        os_command= 'cp m5out/stats.txt stats_injection_rate/stats_{:f}_{:s}.txt'.format(injrate, traffic)
		        os.system(os_command)
                        #print("Threshold local threshold: ")
                        #print(local_threshold)
                    
            # os_command= 'cp m5out/stats.txt stats_injection_rate/stats_{:f}.txt'.format(injrate)
            # os.system(os_command)
            print('final local threshold : ', local_threshold)
            if float(local_threshold) == 100000.0:
                global_threshold = global_threshold
            else:
                global_threshold = local_threshold        
    
            print("Training Phase Ended")
            print("Threshold global threshold: ")
            print(global_threshold)
