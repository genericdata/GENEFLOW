#!/bin/bash

## USAGE: sh merge.sh <fcid> <no_demux: "true" or "false"> <alpha>

fcid=$1
no_demux=$2
alpha=$3

echo "fcid: $fcid"
echo "no demux: $no_demux"
echo "alpha: $alpha"

if [ "$no_demux" == "true" ]; then
    folderType="lane"
else
    folderType="sample"
fi

path="${alpha}/${folderType}/default/${fcid}"
out_dir="${alpha}/merged/${fcid}/merged"

# Ensure output directory exists
mkdir -p "${out_dir}"

# Clear the directory if it already exists with data
rm -f "${out_dir}"/*

# Get number of lanes
num_lanes=$(find "${path}" -mindepth 1 -maxdepth 1 -type d | wc -l)

for file in $(find "${path}/1/" -name "*fastq.gz")
do
    if [ "$no_demux" == "true" ]; then
        id=$(basename "${file%%.fastq.gz}" | cut -d '.' -f2) ################### needs to be same as below
        merged_file="${fcid}_${id}.fastq.gz"
    else
        id=$(perl -ne 'print $1 if /_(n0[1-4]_.+)\.fastq\.gz/' <<< "${file}")
        echo "id: $id"
        merged_file="${fcid}_${id}.fastq.gz"
    fi

    # Loop through the number of lanes to get inputs
    for lane in $(seq 1 $num_lanes); do
        input[$lane]=$(find "${path}/${lane}/" -name "*${id}.fastq.gz")
        echo "input[$lane]: ${input[$lane]}"
    done

    # Merge files
    cat "${input[@]}" > "${out_dir}/${merged_file}"

    echo "Created ${merged_file}"
    echo '--'
done
