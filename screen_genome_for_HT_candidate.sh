#!/bin/bash

# currently set with a Chimp query seq, at 94% coverage (since Krishna uses 94% identity):

#Parameters:
# $1 = Animal (Genome Name)
# $2 = StartValue
# $3 = EndValue

#Example Usage: ./screen_genome_for_HT_candidate.sh Chimp 1 24130

#Make new directory for Genome
cd /mnt/Results/HT_Candidates/L1/
mkdir -p $1
#Make sub-directory to store L1/Genome alignment hits
cd $1
mkdir -p Hits

#Align L1 query seqs to all Genome seqs, using LASTZ
# & and wait indicate the use of multiple cores
cd /data01/Genomes/Vertebrates/$1
for ((i=$2; i<=$3; i++));
do
( lastz /data01/Genomes/Vertebrates/$1/seq$i.fa[unmask,multiple] /home/atma/Testing/Chimp/Chimp_chr2A_106355237_106360249.rev.fasta[unmask] --chain --gapped --coverage=94 --ambiguous=n --ambiguous=iupac --format=general-:name2,start2,end2,score,strand2,size2,name1,start1,end1 > /mnt/Results/HT_Candidates/L1/$1/Hits/LASTZ_L1_$1_seq$i ) &
if (( $i % 30 == 0 )); then wait; fi # Limit to 30 concurrent subshells, so that it doesn't open thousands at once
done
wait

cd /mnt/Results/HT_Candidates/L1/$1/Hits

#Remove all files that are empty
find -size  0 -print0 | xargs -0 rm

#Concatenate all hit files into one file
cat LASTZ_L1_$1_seq* > LASTZ_L1_$1_AllSeqs

#Rearrange columns to put files in BED-like form, to be able to use BEDTools 
for ((i=1; i<=$3; i++));
do
( awk '{print $7 "\t" $8 "\t" $9 "\t" $1 "\t" "1" "\t" $5}' LASTZ_L1_$1_seq$i >> BedFormat_L1_$1_seq$i ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

#For each BED-like file, merge nested or overlapping intervals and create fasta file 
for ((i=1; i<=$3; i++));
do
( mergeBed -s -nms -i BedFormat_L1_$1_seq$i > Merged_L1_$1_seq$i.bed ) &
if (( $i % 30 == 0 )); then wait; fi
done 
wait

for ((i=1; i<=$3; i++));
do
( fastaFromBed -fi /data01/Genomes/Vertebrates/$1/seq$i.fa -bed Merged_L1_$1_seq$i.bed -fo L1_$1_seq$i.fasta ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

#Remove all files that are empty
find -size  0 -print0 | xargs -0 rm

#Concatenate all merged-bed files into one file
cat Merged_L1_$1_seq*.bed > Merged_L1_$1_AllSeqs.bed

#Make a new sub-directory in Genome for strand-correcting files 
cd ..
mkdir Strand_Correct
cd Strand_Correct

#Seperate seqs in each file based on strand
for ((i=1; i<=$3; i++));
do
( grep -h -w '+' ../Hits/Merged_L1_$1_seq$i.bed | fastaFromBed -fi /data01/Genomes/Vertebrates/$1/seq$i.fa -bed stdin -fo Plus_strand_L1_$1_seq$i.fasta ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

for ((i=1; i<=$3; i++));
do
( grep -h -w '-' ../Hits/Merged_L1_$1_seq$i.bed | fastaFromBed -fi /data01/Genomes/Vertebrates/$1/seq$i.fa -bed stdin -fo Minus_strand_L1_$1_seq$i.fasta ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

#Reverse complement seqs on minus strand
for ((i=1; i<=$3; i++));
do
( RevComp Minus_strand_L1_$1_seq$i.fasta RevComp_L1_$1_seq$i.fasta ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

#Concatenate reversed-complemented minus seqs with plus seqs
for ((i=1; i<=$3; i++));
do
( cat RevComp_L1_$1_seq$i.fasta Plus_strand_L1_$1_seq$i.fasta > StrandCorrect_L1_$1_seq$i.fasta ) &
if (( $i % 30 == 0 )); then wait; fi
done
wait

#Remove all files that are empty
find -size  0 -print0 | xargs -0 rm

#Concatenate all strand-correct files into one file
cat StrandCorrect_L1_$1_seq*.fasta > AllSeqs_sc_L1_$1.fasta 
cd ..

#Move this file of strand-correct seqs to a new sub-directory
mkdir Clusters
cd Clusters
cp ../Strand_Correct/AllSeqs_sc_L1_$1.fasta .

#Sort the sequences by length
usearch -sortbylength AllSeqs_sc_L1_$1.fasta -output L1_$1_seqs_sorted_HTcand_94.fasta















