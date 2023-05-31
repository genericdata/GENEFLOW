#!/bin/bash

## USAGE: sh do_fastqc.sh </path/to/data> 

# Load required modules
module purge
module load fastqc/0.11.9
module load multiqc/1.9

# Navigate to directory containing data
cd $1

# Loop over all (lane) directories
for dir in ./*
do
    echo "Processing directory: $dir"

    # Run FastQC
    mkdir -p $dir
    fastqc_files=$(ls $1/${dir}/*.fastq.gz)
    fastqc $fastqc_files -t 80 -o $dir
    
    # Analyze undetermined read barcodes
    echo "Analyzing undetermined read barcodes in: ${1}/${dir}"
    python3 count_barcode_frequency.py $(ls $dir/*_n01_undetermined.fastq.gz) > $dir/undetermined_barcodes.txt
    
    # Create and populate undetermined_barcodes_mqc.txt file
    undetermined_barcodes_mqc=$dir/undetermined_barcodes_mqc.txt
    {
        echo "# id: 'undetermined_barcodes'"
        echo "# plot_type: 'table'"
        echo "# section_name: '	Barcodes of Undetermined Reads'"
        echo "# description: \"<br />We have determined the barcodes of your undetermined reads. Here are the top 20 barcodes. The full list is available <a href='undetermined_barcodes.txt'>here</a>. <b>If your libraries are dual indexed, the two indices are concatenated.</b>\""
        echo "Barcode Sequence(s)	Count	Frequency (%)"
        head -20 $output_dir/undetermined_barcodes.txt
    } > $undetermined_barcodes_mqc

    echo "Finished analyzing undetermined read barcodes"

    # Run MultiQC
    echo "Starting MultiQC"
    cd $dir
    multiqc -f -c /home/gencore/SCRIPTS/hpcpipe/mqc_config.yaml $output_dir
    echo "Finished MultiQC"

    # Return to initial directory
    cd $1
done
