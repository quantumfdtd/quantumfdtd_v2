#!/bin/bash

############################################################################################################
# This script is for multiresolution runs using SLURM
############################################################################################################

############################################################################################################
#SBATCH --job-name=fdtd
#SBATCH --nodes=4
#SBATCH --ntasks 33
#SBATCH --overcommit
############################################################################################################
# number of tasks to use (number of worker cores + 1) 
# should match the SBATCH variable above
tasks=33

# NOTE:  see params.txt file for description of parameter values

# potential
potential=29

# initial condition type
init=1

# reduced quark mass
# bottom
#mass=2.35 
#mass=0.65
mass=1

# phard range
p0=1.0
pf=1.0

# phard steps
psteps=1

# compute deltap
deltap=$(echo "scale=6; ($pf-$p0)/$psteps" | bc)

# list of xi's
xilist=("0" "0.1" "0.5" "1" "2" "5" "8" "10" "13" "18" "20" "30")
nxi=${#xilist[@]}

# initial lattice spacing
a0=0.1

# initial lattice dimensions
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
echo "==> initial phard = " $p0
echo "==> final phard = " $pf
echo "==> psteps = " $psteps
echo "==> deltap = " $deltap
echo "==> xilist = " ${xilist[@]:0}
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
let nxi-=1
for i in `seq 0 $nxi`;
  do
    xi=${xilist[$i]}
    echo "computing xi="$xi
    # reset initial resolutions
    num=$num0
    a=$a0
    eps=$(echo "scale=6; ($a)^2/8" | bc)
    # begin phard loop
    for j in `seq 0 $psteps`;
      do
        p=$(echo "scale=6; $p0+$j*$deltap" | bc)
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
            mpirun -np $tasks ./mpisolve -NUM $num -A $a -EPS $eps -SAVEWAVEFNCS 1 -INITCONDTYPE $init -POTENTIAL $potential -MASS $mass -T $p -XI $xi 
            done  
            # end resolution loop
        else
          printLine
          echo "==> Computing using previous phard wavefunciton"
          printLine
          # no multiresolution needed in this case, just use the previous phard's wavefunction as the initial condition
          mpirun -np $tasks ./mpisolve -NUM $num -A $a -EPS $eps -SAVEWAVEFNCS 1 -INITCONDTYPE 0 -POTENTIAL $potential -MASS $mass -T $p -XI $xi
        fi
        # append output to dat files 
        cat "data/ground_state.out" >> "data/ground_state.dat"
        cat "data/first_excited_state.out" >> "data/first_excited_state.dat"
        cat "data/second_excited_state.out" >> "data/second_excited_state.dat"
      done
      # end phard loop
    done
    # end lxi loop
 
# some cleanup of unnecessary files
rm data/snapshot/*
rm data/wavefunction*
rm data/*.out
rm mpisolve

# replicate the slurm stdout file for the record and then compress it
cp ../../slurm-${SLURM_JOB_ID}.out out.log
gzip out.log

exit
