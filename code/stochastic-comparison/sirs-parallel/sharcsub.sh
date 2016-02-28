#! /bin/bash

# maximum truncation level and
# number of trials at each trunction level
MAXTRUNC=52
NTRIALS=10

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

# output files directory
rdsdir="$basedir/rds"
if [ ! -d $rdsdir ]; then
	mkdir $rdsdir
fi

# sharcnet output files directory
soutdir="$basedir/sout"
if [ ! -d $soutdir ]; then
	mkdir $soutdir
fi

module unload intel
module load r

# create source files and submit
for TRUNC in $(seq 1 $MAXTRUNC); do
	for TRIAL in $(seq 1 $NTRIALS); do
		filebase="fsim-$TRUNC-$TRIAL"
		srcfile="$srcdir/$filebase.r"
		outfile="$outdir/$filebase.Rout"
		cat fsimthread.r | \
			sed -e "s/TRIAL/$TRIAL/" \
				-e "s/TRUNC/$TRUNC/" > $srcfile
		#R CMD BATCH $srcfile $outfile "--args outdir=$rdsdir"
		sqsub -q serial -o "$soutdir/$filebase.%J.out" --mpp 2.5G -r 5h R CMD BATCH "--args outdir=$rdsdir" $srcfile $outfile
	done
done