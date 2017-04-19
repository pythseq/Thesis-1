#!/bin/bash

#Point to list of sequences to obtain free energy and reference fasta file

export LIST=
export RefSeqs=

if [ ! -d Images ]; then
	mkdir Images
fi

#Generate Seq Files

while read line; do
	samtools faidx $RefSeqs "$line" > $line.fa
	RNAfold -p -d2 --noLP < $line.fasta > $line.calc
        cut -f2,3,4 -d" " $line.calc | sed '/TR/d' | sed '/UA/d' | tr -d '[]' | tr -d '()' | tr -d '{' | cut -f1 -d"d"> $line.MFE
	readarray val < $line.MFE
        echo -e $line '\t' ${val[1]} '\t' ${val[2]} '\t' ${val[3]} > $line.fe
	rm $line.fa
	rm $line.calc
	rm $line.MFE
	mv *.ps Images/

done < $LIST

cat *.fe > FREE_ENERGIES.out
