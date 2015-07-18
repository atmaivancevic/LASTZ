#!/bin/bash

# currently set with the HT candidate queries, at 90% coverage (Krishna uses 94% identity):

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
cd /data01/Genomes/Vertebrates/$1

#use parallel --dryrun to see how it's going to look
parallel lastz {} /home/atma/Testing/HT_candidate_queries.fasta[unmask,multiple] --chain --gapped --coverage=70 --identity=90 --ambiguous=n --ambiguous=iupac --format=general-:name2,start2,end2,score,strand2,size2,name1,start1,end1 '>' /mnt/Results/HT_Candidates/L1/$1/Hits/LASTZ_L1_$1_{/.} ::: /data01/Genomes/Vertebrates/$1/seq*.fa

cd /mnt/Results/HT_Candidates/L1/$1/Hits

#Remove all files that are empty
find -size  0 -print0 | xargs -0 rm

#Concatenate all hit files into one file
cat LASTZ_L1_$1_seq* > LASTZ_L1_$1_AllSeqs

#Rearrange columns to put files in BED-like form, to be able to use BEDTools
# & and wait indicate the use of multiple cores, limit to 10 cores 
for ((i=1; i<=$3; i++));
do
( awk '{print $7 "\t" $8 "\t" $9 "\t" $1 "\t" "1" "\t" $5}' LASTZ_L1_$1_seq$i >> BedFormat_L1_$1_seq$i ) &
if (( $i % 10 == 0 )); then wait; fi
done
wait

#For each BED-like file, merge nested or overlapping intervals and create fasta file 
for ((i=1; i<=$3; i++));
do
( mergeBed -s -nms -i BedFormat_L1_$1_seq$i > Merged_L1_$1_seq$i.bed ) &
if (( $i % 10 == 0 )); then wait; fi
done 
wait

for ((i=1; i<=$3; i++));
do
( fastaFromBed -fi /data01/Genomes/Vertebrates/$1/seq$i.fa -bed Merged_L1_$1_seq$i.bed -fo L1_$1_seq$i.fasta ) &
if (( $i % 10 == 0 )); then wait; fi
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
if (( $i % 10 == 0 )); then wait; fi
done
wait

for ((i=1; i<=$3; i++));
do
( grep -h -w '-' ../Hits/Merged_L1_$1_seq$i.bed | fastaFromBed -fi /data01/Genomes/Vertebrates/$1/seq$i.fa -bed stdin -fo Minus_strand_L1_$1_seq$i.fasta ) &
if (( $i % 10 == 0 )); then wait; fi
done
wait

#Reverse complement seqs on minus strand
for ((i=1; i<=$3; i++));
do
( RevComp Minus_strand_L1_$1_seq$i.fasta RevComp_L1_$1_seq$i.fasta ) &
if (( $i % 10 == 0 )); then wait; fi
done
wait

#Concatenate reversed-complemented minus seqs with plus seqs
for ((i=1; i<=$3; i++));
do
( cat RevComp_L1_$1_seq$i.fasta Plus_strand_L1_$1_seq$i.fasta > StrandCorrect_L1_$1_seq$i.fasta ) &
if (( $i % 10 == 0 )); then wait; fi
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
usearch -sortbylength AllSeqs_sc_L1_$1.fasta -output L1_$1_seqs_sorted_HTcand_90.fasta















