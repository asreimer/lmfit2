#!/bin/bash

# get path to make_lmfit2 binary
CWD=`pwd`
LMFIT2=${CWD:0:-6}/bin/make_lmfit2

# dmapdump command (assuming it's in $PATH!)
DMAPDUMP="dmapdump"

# path to expected.lmfit2
EXPECTED=${CWD}/expected.lmfit2
EXPECTEDDUMP=${CWD}/expected.dmapdump
# path to input rawacf file: test.rawacf
INPUT=${CWD}/test.rawacf
# path to test fitted file: output.lmfit2
OUTPUT=${CWD}/output.lmfit2
OUTPUTDUMP=${CWD}/output.dmapdump

# Now run the code and fit the test.rawacf file
echo -e "\nFitting test.rawacf file. This will take a few seconds..."
$LMFIT2 -new $INPUT > $OUTPUT
echo "... done!"
sleep 1
# Now compare the fitted file to expected file
echo "Now comparing output to expected using dmapdump..."
sleep 1
$DMAPDUMP $EXPECTED > $EXPECTEDDUMP
$DMAPDUMP $OUTPUT > $OUTPUTDUMP

#trim out the origin.time and origin.command lines since
#these are expected to be different
sed -i '/origin.time/d' $EXPECTEDDUMP
sed -i '/origin.command/d' $EXPECTEDDUMP
sed -i '/origin.time/d' $OUTPUTDUMP
sed -i '/origin.command/d' $OUTPUTDUMP

# finally use diff to see if files are different
echo -e "Output of a diff:\n"
diff $EXPECTEDDUMP $OUTPUTDUMP

echo -e "\nIf lines are empty, files are the same!"

sleep 1
echo -e "Cleaning up output files...\n"
rm $EXPECTEDDUMP $OUTPUTDUMP $OUTPUT

exit 0

