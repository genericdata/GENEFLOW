#!/bin/bash

## USAGE: sh do_fastqc.sh </path/to/data> <path/to/mqc_config.yaml>

# Load required modules
module purge
module load fastqc/0.11.9
module load multiqc/1.9

# Loop over all (lane) directories
for dir in ${1}/*
do
    echo "Processing directory: $dir"

    # Make and change into directory for output
    dir_name=$(basename "$dir")
    mkdir -p $dir_name
    cd $dir_name

    # Run FastQC
    fastqc_files=$(ls ${dir}/*.fastq.gz)
    fastqc $fastqc_files -t 80 -o .
    
    # Analyze undetermined read barcodes
    echo "Analyzing undetermined read barcodes in: ${dir}"
    count_barcode_frequency.py $(ls $dir/*_n01_undetermined.fastq.gz) > undetermined_barcodes.txt
    
    # Create and populate undetermined_barcodes_mqc.txt file
    undetermined_barcodes_mqc=undetermined_barcodes_mqc.txt
    {
        echo "# id: 'undetermined_barcodes'"
        echo "# plot_type: 'table'"
        echo "# section_name: '	Barcodes of Undetermined Reads'"
        echo "# description: \"<br />We have determined the barcodes of your undetermined reads. Here are the top 20 barcodes. The full list is available <a href='undetermined_barcodes.txt'>here</a>. <b>If your libraries are dual indexed, the two indices are concatenated.</b>\""
        echo "Barcode Sequence(s)	Count	Frequency (%)"
        head -20 undetermined_barcodes.txt
    } > $undetermined_barcodes_mqc

    echo "Finished analyzing undetermined read barcodes"

    # Run MultiQC
    echo "Starting MultiQC"
    mqc_config=$2
    multiqc -f -c $mqc_config .
    echo "Finished MultiQC"

    # Change back to parent directory
    cd ..

done
