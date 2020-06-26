if [ -z $1 ]
then
        echo "usage: processRuns.sh <datekey>"
	exit
fi

# Step in T for fixed XI

XI=0.0
T0=1.0
TF=3.0
TSTEPS=10
TSTEP=$(echo "scale=2; ($TF-$T0)/$TSTEPS" | bc)

DATE=$1
T=$T0
COUNTER=0

OUTPUTFILENAME="sfdtd_bindingenergies_${DATE}_XI_${XI}.dat"

echo "Collecting data in ${OUTPUTFILENAME}"

rm -f $OUTPUTFILENAME

while [ $COUNTER -lt $((TSTEPS+1)) ]; do
  # create directories and configure
  DIRNAME="sfdtd_${DATE}_XI_${XI}_T_${T}"
  cd $DIRNAME
  GSTATE=`awk '/Ground State Binding Energy/' mpisolve.out`
  ESTATE=`awk '/1st excited state binding energy/' mpisolve.out`
  GSTATE=${GSTATE#*:}
  ESTATE=${ESTATE#*:}
  cd ..
  { #output 
    echo "${T}  ${GSTATE}  ${ESTATE}" 
    #
  } >> "${OUTPUTFILENAME}"
  # increment 
  let COUNTER=$COUNTER+1
  T=$(echo "scale=2; $T0 + $COUNTER*$TSTEP" | bc)
done
