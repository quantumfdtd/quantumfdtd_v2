#!/bin/bash

#BE_XI_DEMO  Script to get the binding energy for a specified list of the anisotropy parameter xi.
#
# Synatx:
#	bash be_xi_demo.sh
#
# Out:
#	be_xi{n}.txt - file of the binding energies for xi(n) run
#	xilist.txt   - file of xi parameters to be run.
#	               Format: xi(0) xi(1) ... xi(n)
#
# Description:
# 	Given a specified list of xi parmeters the script sets up a directory
# 	labeled xi with "-" seperating xi(0) ... xi(n) values. PointRun2.py
#	and ExtractBE.py are run for all xi values. The output along with
# 	a file listing the xi values for the runs are dumped into the 
# 	directory.

# Anisotrpy parameter: xi=("xi_0" ... "xi_n")
#xi=("0.0001" "0.1" "0.5" "1" "2" "5" "8" "10" "13" "18" "20" "30")
xi=("0.0001" "1.0" )

# Make directory to store BE files with various xi params.
xidir="xi"
for (( i=0; i<${#xi[@]}; i++))
do
  xidir=$xidir"-${xi[$i]}"
done

# Remove directory if it exits.
if [ -d $xidir ];
then
  rm -rf $xidir
fi
mkdir $xidir

# Output the xi parameters to a file. Usefull for beint.nb.
# xi values are seperated by a single space.
echo ${xi[@]:0} >> $xidir/xilist.txt

# Run the main code.
for (( i=0; i<${#xi[@]}; i++ ));
do
  echo "Running xi = "${xi[$i]}" ..."

  # Clean files.
  ./cleandatafiles.sh
  
  ./PointRun2.py "xi="${xi[$i]}
  ./ExtractBE.py
 
  # Move BE file to xi directory.
  cp be.txt $xidir/be_xi"${xi[$i]}.txt"
done
