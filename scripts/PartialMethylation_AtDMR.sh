#!/bin/bash

echo -e "Arguments are:\n"
echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
echo $7

tabix="$7"
datamash="$6"
if [[ ! -z "$1" && ! -z "$2" && -n "$3" && -n "$4" && ! -z "$5" && ! -z "$6" && ! -z "$7" ]]
then
	echo -e "file_name\tinterval\tNumCG_InFile\tPartialMethylatedCGs\tRetioCGsWithPartialMethylation_DMR\tMeanMethylationAtDMR" > $5 
	for i in `cat $1 | awk '{print $1":"$2"-"$3}'`
	do  
		name=$(echo $2 | rev | cut -d'/' -f1 | rev)
		all=$($tabix $2 $i | wc -l)
		partial=$($tabix $2 $i | awk -v var=$3 '$4 > var' | awk -v var=$4 '$4 < var' | wc -l)
		average=$($tabix $2 $i | $datamash mean 4)
		if [ $all -eq 0 ]
		then
			echo -e "$name\t$i\tNoCG\tNA\tNA\tNA" >> $5 
		else
			echo -e "$name\t$i\t$all\t$partial\t$average" | awk '{print $1,$2,$3,$4,$4/$3,$5}' OFS='\t' >> $5
		fi
	done
else
	echo -e "Please insert positional arguments:\n
	First argument: Path to the DMR file.
	Second argument: Path to bgziped and indexed methylation file. format chr\tstart\tend\tMethylationFrequency
	Third argument: Partial methylation threshold lowerbound (e.g 0.35).
	Forth argument: Partial methylation threshold upperbound (e.g 0.65).
	Fifth argument: Path to the output file.
	Sixth argument: Path to the executable datamash file.
	Seventh argument: Path to the executable tabix file." 
fi
