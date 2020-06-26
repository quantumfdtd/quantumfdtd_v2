if [ -z $1 ]
then
        echo "usage: processRunsVSxi.sh <datekey>"
	exit
fi

# Step in XI for fixed T

T=1.0
XI0=1.0
XIF=3.0
XISTEPS=10
XISTEP=$(echo "scale=2; ($XIF-$XI0)/$XISTEPS" | bc)

DATE=$1
XI=$XI0
COUNTER=0

OUTPUTFILENAME="sfdtd_bindingenergies_${DATE}_T_${T}.dat"

echo "Collecting data in ${OUTPUTFILENAME}"

rm -f $OUTPUTFILENAME

while [ $COUNTER -lt $((XISTEPS+1)) ]; do
  # create directories and configure
  DIRNAME="sfdtd_${DATE}_XI_${XI}_T_${T}"
  cd $DIRNAME
  GSTATE=`awk '/Ground State Binding Energy/' mpisolve.out`
  ESTATE=`awk '/1st excited state binding energy/' mpisolve.out`
  GSTATE=${GSTATE#*:}
  ESTATE=${ESTATE#*:}
  cd ..
  { #output 
    echo "${XI}  ${GSTATE}  ${ESTATE}" 
    #
  } >> "${OUTPUTFILENAME}"
  # increment 
  let COUNTER=$COUNTER+1
  XI=$(echo "scale=2; $XI0 + $COUNTER*$XISTEP" | bc)
done
