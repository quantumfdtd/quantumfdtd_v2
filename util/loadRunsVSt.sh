# Step in T for fixed XI

XI=0.0
T0=1.0
TF=3.0
TSTEPS=10
TSTEP=$(echo "scale=2; ($TF-$T0)/$TSTEPS" | bc)

T=$T0
COUNTER=0
DATE=$(date -u '+%F%H%M%S')
echo $DATE

while [ $COUNTER -lt $((TSTEPS+1)) ]; do
  # create directories and configure
  DIRNAME="sfdtd_${DATE}_XI_${XI}_T_${T}"
  cp -rf ./template "$DIRNAME"
  cd $DIRNAME
  SEDRULE="s/sfdtd/${DIRNAME}/g"
  sed $SEDRULE runner.pbs >/tmp/$$ && mv /tmp/$$ runner.pbs
  SEDRULE="s/temperature/${T}/g"
  sed $SEDRULE runner.pbs >/tmp/$$ && mv /tmp/$$ runner.pbs
  SEDRULE="s/anisoparameter/${XI}/g"
  sed $SEDRULE runner.pbs >/tmp/$$ && mv /tmp/$$ runner.pbs
  qsub -q myri runner.pbs
  cd ..
  # increment 
  let COUNTER=$COUNTER+1
  T=$(echo "scale=2; $T0 + $COUNTER*$TSTEP" | bc)
done
