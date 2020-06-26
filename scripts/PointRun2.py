#!/usr/bin/python

# POINTRUN2  
#
#  Syntax:
#    PointRun2.py xi=5
#
#  In:
#    XI - Anisotropy parameter   (optional, default 5)
#
#  Description:
#    Extension of PointRun.py to take imports from the system.

import os
import popen2
import sys

# ==================================================================
# === Begin Variable Declaration ===================================
# ==================================================================

# Potential to Use
potential = 21

# Mass to use (Charmonium = 0.65, Bottomonium = 2.35)
#mass = 0.65
mass = 2.35

# Anisotropy factor xi
xi = 5

# Temperature Range
tStart = 1.0
tFinish = 4.0

# Lattice Size Factor Maximum (jMax < 4!)
jMax = 2

# Number of Nodes (n = [Servers] * 8 + 1)
n = 33

#
# Lattice Spacing
#
# j/psi and upsilon
#aStart = 0.4 
#
# xib
aStart = 0.6 

# Step Size and Step Count
nSteps = 30
delta = (tFinish - tStart) / nSteps

# Import from system.
for arg in sys.argv:
	var = arg.partition("=")
	if var[0] ==  "xi":
		xi = float(var[2])		

# ==================================================================
# === Begin Method =================================================
# ==================================================================
os.system('echo Initiating...')

# Run the method for each point between the temperatures
for i in range (0, nSteps + 1):
	# Get the temperature at the current point
	temp = tStart + i * delta
	print "Running temperature %f" % temp
	
	# Run the points starting with a rough lattice and get finer, using the
	# rough lattice as a base for each successive finer lattice
	if i==0: 
		for j in range (0, jMax + 1):
			# Calculate the number of lattice spaces
			num = 2 * (n - 1) * (2 ** j)
		
			# Determine whether the wavefunctions should be saved
			save = 1 # always save!
		
			# Determine if the wavefunctions should be loaded or generated randomly		
			init = 1
			if j > 0:
				init = 0
			
			# Calculate the needed lattice spacing
			a = aStart / (2 ** j)
			eps = (a * a) / 8
		
			print "Current Box Size: %f; Nodes: %f" % (num, n)
		
			# Generate the name for the output file
			boxNum = "%s_%s" % (i,num)
			if num == (2 * (n - 1) * (2 ** jMax)):
				boxNum = "%s" % i
			# Run the code
			os.system("mpirun --hostfile ~/my-hostfile.txt -np %s ./mpisolve -T %s -MASS %s -XI %s -POTENTIAL %s  -INITCONDTYPE %s -SAVEWAVEFNCS %s -NUM %s -EPS %s -A %s > zrun%s.log 2> zrunErr.log" % (n,temp,mass,xi,potential,init,save,num,eps,a,boxNum))
	else:
			num = 2 * (n-1) * (2 ** jMax)	
			print "Current Box Size: %f; Clusters: %f" % (num, n)
			save = 1
			init = 0
			a = aStart / (2 ** jMax)
			eps = (a * a)/8
			boxNum = "%s" % i
			# Run the code
			os.system("mpirun --hostfile ~/my-hostfile.txt -np %s ./mpisolve -T %s -MASS %s -XI %s -POTENTIAL %s  -INITCONDTYPE %s -SAVEWAVEFNCS %s -NUM %s -EPS %s -A %s > zrun%s.log 2> zrunErr.log" % (n,temp,mass,xi,potential,init,save,num,eps,a,boxNum))
		
