# Step in T for fixed XI

T=1.0
XI0=1.0
XIF=3.0
XISTEPS=10
XISTEP=$(echo "scale=2; ($XIF-$XI0)/$XISTEPS" | bc)

XI=$XI0
COUNTER=0
DATE=$(date -u '+%F%H%M%S')
echo $DATE

while [ $COUNTER -lt $((XISTEPS+1)) ]; do
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
  XI=$(echo "scale=2; $XI0 + $COUNTER*$XISTEP" | bc)
done
