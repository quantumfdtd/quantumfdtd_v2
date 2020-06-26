#!/bin/bash

############################################################################################################
# This script is for multiresolution runs using SLURM
############################################################################################################

############################################################################################################
#SBATCH --job-name=fdtd
#SBATCH --nodes=4
#SBATCH --ntasks 33
#SBATCH --overcommit
#SBATCH --exclusive
############################################################################################################
# number of tasks to use (number of worker cores + 1) 
# should match the SBATCH variable above
tasks=33

# NOTE:  see params.txt file for description of parameter values

# potential
# charm
# potential=32
# shifted charm
potential=31
# bottom
# potential=33

# initial condition type
init=3

# reduced quark mass
# bottom
#mass=2.35 
# charm
mass=0.645

# eB range
eB0=0.0
eBf=0.3

# eB steps
eBsteps=6

# compute deltap
deltaeB=$(echo "scale=6; ($eBf-$eB0)/$eBsteps" | bc)

# OVERRIDE!
# eBsteps=0

# list of kx's
kxlist=("0" "0.1" "0.25" "0.5" "0.75" "1" "2" "3" "4" "5")
#kxlist=("0" "0.1" "0.25" "0.5" "0.75"  "1" "2" "5" "7" "10")
#kxlist=("7")
#kxlist=("10")
nkx=${#kxlist[@]}

# initial lattice spacing
a0=0.4

# initial lattice dimensions
# upsilon
#num0=128
# j/psi
num0=128

# number of multiresolution steps to make
ne=2

############################################################################################################
############################################################################################################

function printLine {
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}

printLine
echo "Initializing variables"
printLine
echo "==> potential = " $potential
echo "==> reduced quark mass = " $mass
echo "==> initial eB = " $eB0
echo "==> final eB = " $eBf
echo "==> eBsteps = " $eBsteps
echo "==> deltaeB = " $deltaeB
echo "==> kxlist = " ${kxlist[@]:0}
echo "==> initial lattice spacing = " $a0
echo "==> initial lattice dimension = " $num0
printLine

# check for existence of 'runs' directory; create it if it doesn't exist
if [ -d 'runs' ]; then
	echo "==> 'runs' directory exists"
else 
	echo "==> 'runs' directory does not exist, creating it"
	mkdir runs
fi

# check for existence of run specfic directory
id=`eval date +%Y%m%d`"_"$SLURM_JOB_ID
dir="runs/run_"$id
if [ -d $dir ]; then
	echo "==> run directory "$dir" directory exists ... exiting"
	exit
else 
	echo "==> creating run directory "$dir
	mkdir $dir
fi

# copy necessary files and create output directories
echo "==> creating run directory structure"

mkdir $dir"/data"
mkdir $dir"/data/snapshot"
cp mpisolve $dir
cp params.txt $dir
cd $dir

printLine
echo "==> slurm job id "$SLURM_JOB_ID
echo "==> slurm node list "$SLURM_JOB_NODELIST
echo "==> slurm tasks per node "$SLURM_TASKS_PER_NODE
printLine
echo "==> Starting computation..."
printLine

# now, finally, execute the mpi jobs
let nkx-=1
for i in `seq 0 $nkx`;
  do
    kx=${kxlist[$i]}
    echo "computing kx="$kx
    # reset initial resolutions
    num=$num0
    a=$a0
    eps=$(echo "scale=6; ($a)^2/8" | bc)
    # begin dB loop
    for j in `seq 0 $eBsteps`;
      do
        eB=$(echo "scale=6; $eB0+$j*$deltaeB" | bc)
        if [ $j -eq 0 ]; then
          # begin multiresolution loop
          for k in `seq 1 $ne`;
            do
              if [ $k -gt 1 ]; then 
                init=0
                a=$(echo "scale=6; $a/2" | bc)
                num=$(echo "scale=0; $num*2" | bc)
                eps=$(echo "scale=6; $eps/4" | bc)
              fi
            printLine
            echo "==> Computing at resolution "$k", num="$num", a="$a", eps="$eps
            printLine
            mpirun -np $tasks ./mpisolve -NUM $num -A $a -EPS $eps -SAVEWAVEFNCS 1 -INITCONDTYPE $init -POTENTIAL $potential -MASS $mass -eB $eB -Kx $kx 
            done  
            # end resolution loop
        else
          printLine
          echo "==> Computing using previous eB wavefunciton"
          printLine
          # no multiresolution needed in this case, just use the previous wavefunction as the initial condition
          mpirun -np $tasks ./mpisolve -NUM $num -A $a -EPS $eps -SAVEWAVEFNCS 1 -INITCONDTYPE 0 -POTENTIAL $potential -MASS $mass -eB $eB -Kx $kx 
        fi
        # append output to dat files 
        cat "data/ground_state.out" >> "data/ground_state.dat"
        cat "data/first_excited_state.out" >> "data/first_excited_state.dat"
        cat "data/second_excited_state.out" >> "data/second_excited_state.dat"
      done
      # end eB loop
    done
    # end kx loop
 
# some cleanup of unnecessary files
rm data/snapshot/*
rm data/wavefunction*
rm data/*.out
rm mpisolve

# replicate the slurm stdout file for the record and then compress it
cp ../../slurm-${SLURM_JOB_ID}.out out.log
gzip out.log

exit
