#! /bin/bash

# output files directory
outdir="out"
if [ ! -d $outdir ]; then
	mkdir $outdir
fi

# output files directory
logdir="log"
if [ ! -d $logdir ]; then
	mkdir $logdir
fi

for datafile in $(ls data/data*.txt); do

	threadnum=$(echo $datafile | sed -e s/[^0-9]//g)

	neinumfile="data/neinum-$threadnum.txt"
	neibmatfile="data/neibmat-$threadnum.txt"
	outfile="$outdir/cuIF2states-$threadnum.dat"
	logfile="$logdir/cuIF2-$threadnum.log"

	./cuIF2 $datafile $neinumfile $neibmatfile $outfile > $logfile

done
