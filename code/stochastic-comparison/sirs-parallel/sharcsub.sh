#! /bin/bash

# maximum truncation level and
# number of trials at each trunction level
NTRIALS=8

basedir="/work/barrowdd"

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
if [ ! -d $varimagesdir ]; then
	mkdir $varimagesdir
fi

# sharcnet output files directory
soutdir="$basedir/sout"
if [ ! -d $soutdir ]; then
	mkdir $soutdir
fi

module unload intel
module load r

# create source files and submit
for TRIAL in $(seq 1 $NTRIALS); do
	filebase="fsim-$TRIAL"
	srcfile="$srcdir/$filebase.r"
	outfile="$outdir/$filebase.Rout"
	cat fsimthread.r | sed -e "s/TRIAL/$TRIAL/" > $srcfile
	sqsub -q serial -o "$soutdir/$filebase.%J.out" --mpp 2.5G -r 20h Rscript $srcfile $vardir
done
