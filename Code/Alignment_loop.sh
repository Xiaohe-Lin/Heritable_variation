#!/bin/bash
# Make sure you have already installed trimmomatic and star, and prepared all the index and files needed including star index，genome and gtf files.
# Function to process a pair of files using trimmomatic and STAR
process_files() {
    file1="$1"
    file2="$2"
    Root="$3"
    result="/transcript_analyze/STAR/result/"

    echo "#####---Remove adapter and low quality reads---#####"
    trimmomatic PE -threads 15 $file1 $file2 $result$Root'_1.paired.fq.gz' $result$Root'_1.unpaired.fq.gz' $result$Root'_2.paired.fq.gz' $result$Root'_2.unpaired.fq.gz' "ILLUMINACLIP:/home/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    rm $result$Root'_1.unpaired.fq.gz'
    rm $result$Root'_2.unpaired.fq.gz'

    echo "#####---Mapping to genome---#####"
    STAR --runThreadN 15 --runMode alignReads --genomeDir /transcript_analyze/STAR/index --readFilesCommand zcat --sjdbGTFfile TAIR10.gtf --sjdbOverhang 149 --readFilesIn "$result$Root'_1.paired.fq.gz'" "$result$Root'_2.paired.fq.gz'" --outFileNamePrefix "$result$Root" --outSAMtype BAM SortedByCoordinate
}

# Set the directory containing the fastq files
fastq_dir="/path/to/folder/with/fastq/files"

# Iterate over all the fastq.gz files in the directory
for file1 in "$fastq_dir"/*_1.fq.gz; do
    # Get the corresponding file2 based on file1
    file2="${file1/_1.fq.gz/_2.fq.gz}"
    
    # Extract the Root name from file1 (assuming it follows a certain naming pattern)
    Root=$(basename "$file1" | cut -d '_' -f 1)
    
    echo "####---Processing $file1 and $file2...---#####"
    
    # Call the function to process the pair of files
    process_files "$file1" "$file2" "$Root"
done
