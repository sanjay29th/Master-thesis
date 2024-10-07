#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 input_fasta_file output_fasta_file"
    exit 1
}

# Check if correct number of arguments is given
if [ $# -ne 2 ]; then
    usage
fi

input_file=$1
output_file=$2

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file not found!"
    exit 1
fi

# Initialize variables
contig_number=0
seq_length=0
sequence=""

# Create or clear the output file
> "$output_file"

# Function to write the current contig's header and sequence to the output file
write_contig() {
    if [[ -n "$sequence" ]]; then
        echo ">NODE_${contig_number}_length_${seq_length}_cov_0" >> "$output_file"
        echo "$sequence" >> "$output_file"
    fi
}

# Read the input FASTA file and process it
while IFS= read -r line; do
    if [[ $line =~ ^\> ]]; then
        # Write the previous contig to the output file
        write_contig

        # Start a new contig
        contig_number=$((contig_number + 1))
        seq_length=0
        sequence=""
    else
        # Append the line to the current sequence
        sequence+="$line"
        seq_length=$((${#line} + seq_length))
    fi
done < "$input_file"

# Write the last contig to the output file
write_contig

echo "Headers rewritten and output saved to $output_file"
