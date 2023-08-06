#!/bin/bash
#---------------------------------#
# 1.Quality control and alignment #
#---------------------------------#

function qc() {
	trimmomatic PE -phred33 $1"_1.fq.gz" $1"_2.fq.gz" $1"_1_paired.fq.gz" $1"_1_unpaired.fq.gz" $1"_2_paired.fq.gz" $1"_2_unpaired.fq.gz" "ILLUMINACLIP:/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
	fastqc -t 8 -f fastq $1"_1.fq.gz" $1"_2.fq.gz"
	fastqc -t 8 -f fastq $1"_1_paired.fq.gz" $1"_2_paired.fq.gz"
	bwameth.py index $2".fa"
	bwameth.py --reference $2".fa" $1"_1_paired.fq.gz" $1"_2_paired.fq.gz" -t 12 | samtools view -b -> sample.bwameth.bam
    sample_bam=sample.bwameth.bam
}

#---------------------------------------------------#
# 2.Filering alignment results (mapping quality>30) #
#---------------------------------------------------#
function filter(){
	samtools view -@ 24 -b -q 30 $1 > sample.q30.bam
	sample_q30_bam=sample.q30.bam
}

#-----------------------------#
# 3.Sorting alignment results #
#-----------------------------#
function sort () {
	samtools sort -o sample.sorted.bam $1
    sorted_bam=sample.sorted.bam
}


#----------------------------------#
# 4.Duplicates marked and removerd #
#----------------------------------#
function markdup () {
	picard MarkDuplicates \
            -I $1 \
            -O sample.markdup.bam \
            -M sample.markdup.txt 
    samtools view -b -F 1024 -h sample.markdup.bam > sample.rmdup.bam
    markdup_bam=sample.markdup.bam
    rmdup_bam=sample.rmdup.bam
}


#------------------------------------------------#
# 5.extracting and filtering methylation metrics #
#------------------------------------------------#
function extract(){
	samtools index $1
	MethylDackel extract $2".fa" $1 -p 20 -q 40 -o sample --cytosine_report --CHG --CHH  --OT 0,0,0,135 --OB 0,0,5,0
}

#---------------------------#
# 6.package above functions #
#---------------------------#

function bwameth () {
	qc $1 $2 
	rm $1"_1_unpaired.fq.gz" $1"_2_unpaired.fq.gz"
	filter $sample_bam
	sort $sample_q30_bam
	markdup $sorted_bam
	extract $rmdup_bam $2
}

bwameth $1 $2
