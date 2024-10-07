#!/bin/bash
set -euo pipefail

# Initialize Conda (make sure the path is correct)
source ~/miniconda3/etc/profile.d/conda.sh

# Function to display help message
show_help() {
    echo "Usage: $0 -i input_fasta_file"
    echo "  -i input_fasta_file   Specify the input genome fasta sequence file."
}

# Parse command-line arguments
input_fasta=""
while getopts ":i:h" opt; do
    case ${opt} in
        i )
            input_fasta="${OPTARG}"
            ;;
        h )
            show_help
            exit 0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            show_help
            exit 1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" >&2
            show_help
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

# Check if the input fasta file is provided
if [ -z "$input_fasta" ]; then
    echo "Error: Input fasta file is required." >&2
    show_help
    exit 1
fi

# Define the base name for the output directories and files
base_name=$(basename "$input_fasta" .fasta)

# Ensure output directory exists
output_dir="${base_name}_output"
mkdir -p "$output_dir"
echo "Created output directory: $output_dir"

# Step 1: Run VirSorter with the input fasta file
echo "Running VirSorter..."
conda activate virsorter2_env || { echo "Error: Conda environment 'virsorter2_env' activation failed." >&2; exit 1; }
virsorter run -w "${output_dir}/virsorter2_output/" -i "$input_fasta" --include-groups "dsDNAphage,ssDNA,RNA" -j 8 -d /mnt/colinb/databases/virsorter2/
conda deactivate
if [ $? -ne 0 ]; then
    echo "Error: VirSorter failed." >&2
    exit 1
fi
echo "VirSorter completed successfully."

# Step 2: Run CheckV
echo "Running CheckV..."
conda activate checkv_env || { echo "Error: Conda environment 'checkv_env' activation failed." >&2; exit 1; }
checkv end_to_end "${output_dir}/virsorter2_output/final-viral-combined.fa" "${output_dir}/checkv_output" -d /mnt/sanjay/db/checkv-db-v1.5 -t 8
conda deactivate
if [ $? -ne 0 ]; then
    echo "Error: CheckV failed." >&2
    exit 1
fi
echo "CheckV completed successfully."

# Step 3: Annotate prophage genomes with Pharokka
echo "Annotating prophage genomes with Pharokka..."
conda activate pharokka_env || { echo "Error: Conda environment 'pharokka_env' activation failed." >&2; exit 1; }
pharokka.py -i "${output_dir}/checkv_output/proviruses.fna" -o "${output_dir}/pharokka_output" -d /mnt/colinb/databases/pharokka/ -t 8
conda deactivate
if [ $? -ne 0 ]; then
    echo "Error: Pharokka annotation failed." >&2; exit 1
fi
echo "Pharokka annotation completed successfully."

echo "Pipeline steps completed successfully. Proceeding to data transformation..."

# Step 4: Extract details of interest to obtain required tsv files and prophage sequences

# Extract the first column from final-viral-score.tsv
echo "Extracting details from final-viral-score.tsv..."
cut -f1 "${output_dir}/virsorter2_output/final-viral-score.tsv" > "${output_dir}/virsorter2_output/temp-column.txt"

# Combine the extracted column with the rest of final-viral-boundary.tsv
echo "Combining extracted column with final-viral-boundary.tsv..."

awk 'NR==1 ||  $1=$29' OFS="\t" "${output_dir}/virsorter2_output/final-viral-boundary.tsv" > "${output_dir}/virsorter2_output/final-viral-boundary-new.tsv"

# Clean up temporary file
rm "${output_dir}/virsorter2_output/temp-column.txt"

# Merge final-viral-score.tsv and final-viral-boundary-new.tsv to obtain boundary details
echo "Merging final-viral-score.tsv and final-viral-boundary-new.tsv..."

# Define paths to input files
score_file="${output_dir}/virsorter2_output/final-viral-score.tsv"
boundary_file="${output_dir}/virsorter2_output/final-viral-boundary-new.tsv"
output_file="${output_dir}/virsorter2_output/merged_final_viral_data.tsv"

# Echo the header into a temporary file
echo -e "seqname\tdsDNAphage\tssDNA\tRNA\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular\ttrim_bp_start\ttrim_bp_end" > "${output_file}.tmp"

# Perform the join operation and append the results to the temporary file
join -t $'\t' \
    <(sort -k1,1 "$score_file") \
    <(awk -F'\t' 'NR>1 {print $1 "\t" $4 "\t" $5}' "$boundary_file" | sort -k1,1) \
    >> "${output_file}.tmp"

# Check if join was successful and if the output file has been created
if [ $? -eq 0 ] && [ -s "${output_file}.tmp" ]; then
    mv "${output_file}.tmp" "$output_file"
    echo "Merging of final-viral-score.tsv and final-viral-boundary-new.tsv completed successfully."
else
    echo "Error: Join operation failed or output file is empty."
    rm -f "${output_file}.tmp"
    exit 1
fi

# Extract header details from checkv proviruses.fna
echo "Extracting header details from checkv proviruses.fna..."
grep '^>' "${output_dir}/checkv_output/proviruses.fna" | sed 's/^>//; s/[ \t\/-]/\t/g' | { echo -e "seqname\tcheckv_left\tcheckv_right\tcheckv_length"; cat -; } > "${output_dir}/checkv_output/temp" && mv "${output_dir}/checkv_output/temp" "${output_dir}/checkv_output/checkv_headers.tsv"

# Merge merged_final_viral_data.tsv and checkv_headers.tsv based on the first column
echo "Merging merged_final_viral_data.tsv and checkv_headers.tsv..."

# Define file paths
viral_data="${output_dir}/virsorter2_output/merged_final_viral_data.tsv"
checkv_data="${output_dir}/checkv_output/checkv_headers.tsv"
checkv_data_modified="${output_dir}/checkv_output/checkv_headers_modified.tsv"
output_file="${output_dir}/final_merged_data.tsv"

# Extract headers
viral_header=$(head -n 1 "$viral_data")
checkv_header=$(head -n 1 "$checkv_data" | cut -f2-)

# Modify checkv_data to remove suffixes
awk -F'\t' '{if(NR>1) {sub(/_[0-9]+$/, "", $1);} print}' OFS='\t' "$checkv_data" > "$checkv_data_modified"

# Perform join operation, excluding headers
join -t $'\t' \
    <(tail -n +2 "$viral_data" | sort -k1,1) \
    <(tail -n +2 "$checkv_data_modified" | sort -k1,1) \
    > "$output_file.tmp"

# Combine headers and prepend to output
echo -e "$viral_header\t$checkv_header" > "$output_file"
cat "$output_file.tmp" >> "$output_file"

# Clean up temporary file
rm "$output_file.tmp"
echo "Merging of final_merged_data.tsv completed successfully."

# Determine prophage boundaries on host genomes
echo "Determining prophage boundaries on host genomes..."
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "prophage_left", "prophage_right"} NR>1 {prophage_left = $11 + $13; prophage_right = $11 + $14; print $0, prophage_left, prophage_right}' "${output_dir}/final_merged_data.tsv" > "${output_dir}/script_output.tsv"

# Generate manual_prophage_coordinates file for propagate
echo "Generating manual_prophage_coordinates file..."

# Input TSV file name
input_file="${output_dir}/script_output.tsv"

# Output TSV file name
output_file="${output_dir}/manual_prophage_coordinates.tsv"

# Initialize counter for prophage numbering
prophage_counter=1

# Process the file with awk
awk 'BEGIN { FS="\t"; OFS="\t"; prophage_counter=0; }
    {
        # Duplicate seqname to scaffold
        scaffold = $1;
        # Remove ||*_partial and ||*_full from scaffold column
        gsub(/\|\|[^|]*_partial|\|\|full/, "", scaffold);
        # Increment prophage counter and format fragment column
        if (NR > 1) {
            fragment = "pharokka_" ++prophage_counter;
        } else {
            fragment = "fragment";  # header line
        }

        # Print the modified line with new headers
        if (NR == 1) {
            print "scaffold", "fragment", "start", "stop";
        } else {
            print scaffold, fragment, $16, $17;
        }
    }
' "$input_file" > "$output_file"

echo "Transformation complete. Output saved to $output_file"
