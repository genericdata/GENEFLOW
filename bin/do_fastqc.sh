#!/bin/bash

## USAGE: sh do_fastqc.sh <path_to_data>
## EXAMPLE: sh do_fastqc.sh /scratch/mk5636/alpha/lane/HH2KLM7XX/1
## EXAMPLE: sh do_fastqc.sh /scratch/mk5636/alpha/sample/HH2KLM7XX/1
## EXAMPLE: sh do_fastqc.sh /scratch/mk5636/alpha/merged/HWDFCAKXY/merged

# Load required modules
module purge
module load fastqc/0.11.9

echo "Starting FastQC"

# Get directory to analyze
dir=$1

# Run FastQC
fastqc_files=$(ls ${dir}/*.fastq.gz)
fastqc $fastqc_files -t 80 -o .

# Analyze undetermined read barcodes
file=$(ls $dir/*_n01_undetermined.fastq.gz 2> /dev/null)
if [ -e "$file" ]
then
    echo "Analyzing undetermined read barcodes in: ${dir}"
    count_barcode_frequency.py "$file" > undetermined_barcodes.txt

    # Create and populate undetermined_barcodes_mqc.txt file
    undetermined_barcodes_mqc=undetermined_barcodes_mqc.txt
    {
        echo "# id: 'undetermined_barcodes'"
        echo "# plot_type: 'table'"
        echo "# section_name: '	Barcodes of Undetermined Reads'"
        echo "# description: \"<br />We have determined the barcodes of your undetermined reads. Here are the top 20 barcodes. The full list is available <a href='undetermined_barcodes.txt'>here</a>.\""
        head -21 undetermined_barcodes.txt
    } > $undetermined_barcodes_mqc

    echo "Finished analyzing undetermined read barcodes"
fi

