#!/bin/bash

# This script takes three gfa files: unitig gfa, contig hap1 gfa, and contig hap2 gfa
# With a provided unitig name, the program will find the corresponding contig name that has this unitig
# Vice versa, if a contig name is provided, the program will find all unitigs that contributes to the assembly of this contig
#
# Optional:
# If users want to find where the unitigs and contigs have been aligned to, input PAF alignment files.
# PAF lines containing the alignment blocks will be returned.

# Functions
findctg () {
	reads="$(grep -w "A" $utgFile | grep $unitig | awk '{print $5}')"
	local prev_ctg="ctg"
	local curr_ctg="ctg"
	for i in $reads
	do
		tmp="$(grep "$i" $ctgFile1 $ctgFile2 | awk '{print $2}')"
		if [[ $tmp = "" ]]; then continue; fi
		curr_ctg="$tmp"
		if [[ "$prev_ctg" != "ctg" && "$curr_ctg" != "$prev_ctg" ]]; then
			echo "Caution: Read $i from unitig $unitig is in both $prev_ctg and $curr_ctg."
		fi
		prev_ctg="$curr_ctg"
	done
	contig="$curr_ctg"
}

findutgs () {
	echo "Finding reads in this contig..."
	reads="$(grep -w "A" $ctgFile1 $ctgFile2 | grep $contig | awk '{print $5}')"
	read_cnt="$(echo "$reads"|wc -l)"
	echo "Read count: $read_cnt"
	echo "Finding unitigs with these reads..."
	unitig="$(for i in "$reads"; do grep "$i" $utgFile | awk '{print $2}'; done | uniq)"
	echo "Done."
}

findalns () {
	queries="$1"
	alnFile=$2
	local alnBlocks="$(for i in $queries; do grep "$i" $alnFile | grep tp:A:P; done)"
	echo "$alnBlocks"
}

# Begin running script...
# echo "Unitig GFA file: $1"
# echo "contig hap1 GFA file: $2"
# echo "contig hap2 GFA file: $3"

utgFile="$1"
ctgFile1="$2"
ctgFile2="$3"
read -p "Query unitig name (press enter if proceeding to query contig): " unitig
if [[ "$unitig" = '' ]]
then
	# Enter hit, move to contig name user input
	read -p "Query contig name: " contig
	findutgs
	echo "The unitigs from contig $contig: ${unitig//$'\n'/,}"
else
	# Find the corresponding contig name that has this unitig
	findctg
	echo "The contig with unitig $unitig: $contig"
fi

# Optional: find where the unitigs and contigs have been aligned to and input PAF alignment files
read -p "Keep finding where the unitig(s) and contig(s) have been aligned to? (y/n): " -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
	echo
	read -p "Unitig alignment file: " utgAlnFile
	if [[ $contig == *"h1"* ]]; then
		read -p "Contig alignment file (should be hap1): " ctgAlnFile
	else
		read -p "Contig alignment file (should be hap2): " ctgAlnFile
	fi
	cwd=$(pwd)
	echo "Finding query unitig(s) in $utgAlnFile ..."
	findalns "$unitig" "$utgAlnFile" > query.utg.aln.paf
    echo "Done. Saving output to $cwd/query.utg.aln.paf"

    echo "Finding query contig(s) in $ctgAlnFile ..."
	findalns "$contig" "$ctgAlnFile" > query.ctg.aln.paf
    echo "Done. Saving output to $cwd/query.ctg.aln.paf"
else
	echo
	[[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 
fi