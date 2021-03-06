#! /bin/bash -e

# maximum truncation level and
# number of trials at each trunction level
TRIALMIN=1
TRIALMAX=20

## base directory
basedir="orca"
if [ ! -d  $basedir ]; then
	mkdir $basedir
fi

# source files directory
srcdir="$basedir/src"
if [ ! -d  $srcdir ]; then
	mkdir $srcdir
fi

# output files directory
outdir="$basedir/out"
if [ ! -d $outdir ]; then
	mkdir $outdir
fi

# image files directory
vardir="$basedir/varimages"
if [ ! -d $vardir ]; then
	mkdir $vardir
fi

# sharcnet output files directory
soutdir="$basedir/sout"
if [ ! -d $soutdir ]; then
	mkdir $soutdir
fi

module unload intel
module load r

# create source files and submit
for TRIAL in $(seq $TRIALMIN $TRIALMAX); do
	filebase="fsim-$TRIAL"
	srcfile="$srcdir/$filebase.r"
	outfile="$outdir/$filebase.Rout"
	cat fsimthread.r | sed -e "s/TRIAL/$TRIAL/" > $srcfile
	Rscript $srcfile $vardir > "$soutdir/$filebase.out" 2>&1 &
	echo "Job $TRIAL submitted, sleeping ..."
	sleep 60
done
