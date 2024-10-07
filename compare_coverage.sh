#!/bin/bash

# Usage: ./compare_coverage.sh <forward_reads> <reverse_reads> <genome_fasta> <coordinates_file> <output_file> <threads> <output_dir>

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <forward_reads> <reverse_reads> <genome_fasta> <coordinates_file> <output_file> <threads> <output_dir>"
    exit 1
fi

FORWARD_READS=$1
REVERSE_READS=$2
GENOME_FASTA=$3
COORD_FILE=$4
OUTPUT_FILE=$5
THREADS=$6
OUTPUT_DIR=$7

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define output file names
SAM_FILE="$OUTPUT_DIR/$(basename "$FORWARD_READS" .fastq).sam"
UNSORTED_BAM_FILE="$OUTPUT_DIR/$(basename "$FORWARD_READS" .fastq).unsorted.bam"
SORTED_BAM_FILE="$OUTPUT_DIR/$(basename "$FORWARD_READS" .fastq).sorted.bam"

# Index the genome
echo "Indexing the genome with bowtie2-build..."
bowtie2-build "$GENOME_FASTA" "$OUTPUT_DIR/genome_index"

# Align reads to the genome and create SAM file
echo "Aligning reads with bowtie2..."
echo "SAM output will be saved to: $SAM_FILE"
bowtie2 -x "$OUTPUT_DIR/genome_index" -1 "$FORWARD_READS" -2 "$REVERSE_READS" -S "$SAM_FILE" -p "$THREADS"

# Check if SAM file was generated
if [ ! -s "$SAM_FILE" ]; then
    echo "Error: SAM file was not generated"
    exit 1
fi

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -bS "$SAM_FILE" > "$UNSORTED_BAM_FILE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to convert SAM to BAM"
    exit 1
fi

# Check if BAM file was generated
if [ ! -s "$UNSORTED_BAM_FILE" ]; then
    echo "Error: Unsorted BAM file was not generated"
    exit 1
fi

# Sort BAM file
echo "Sorting BAM file..."
samtools sort "$UNSORTED_BAM_FILE" -o "$SORTED_BAM_FILE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to sort BAM file"
    exit 1
fi

# Check if sorted BAM file was generated
if [ ! -s "$SORTED_BAM_FILE" ]; then
    echo "Error: Sorted BAM file was not generated"
    exit 1
fi

# Index the sorted BAM file
echo "Indexing sorted BAM file..."
samtools index "$SORTED_BAM_FILE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to index BAM file"
    exit 1
fi

# Check if index was generated
if [ ! -s "${SORTED_BAM_FILE}.bai" ]; then
    echo "Error: BAM index file was not generated"
    exit 1
fi

# Optionally, clean up unsorted BAM file if no longer needed
# rm -f "$UNSORTED_BAM_FILE"

echo "Alignment and processing completed successfully!"

# Write the header to the output file
echo -e "scaffold\tfragment\tprophage_mean_coverage\tprophage_median_coverage\tentire_scaffold_mean_coverage\tentire_scaffold_median_coverage\tnon_prophage_mean_coverage\tnon_prophage_median_coverage\tmean_coverage_ratio\tmedian_coverage_ratio" > "$OUTPUT_FILE"

# Function to calculate mean coverage for a region
calculate_mean_coverage() {
    local scaffold=$1
    local start=$2
    local end=$3
    
    samtools depth -r "${scaffold}:${start}-${end}" "$SORTED_BAM_FILE" | \
    awk '{sum += $3} END {if (NR > 0) print sum / NR; else print 0}'
}

# Function to calculate median coverage for a region
calculate_median_coverage() {
    local scaffold=$1
    local start=$2
    local end=$3
    
    samtools depth -r "${scaffold}:${start}-${end}" "$SORTED_BAM_FILE" | \
    awk '{print $3}' | sort -n | awk '{
        count[NR] = $1;
        total = NR;
    } END {
        if (total % 2 == 1) {
            print count[(total + 1) / 2];
        } else {
            print (count[total / 2] + count[total / 2 + 1]) / 2;
        }
    }'
}

# Read the coordinates file line by line, skipping the header
tail -n +2 "$COORD_FILE" | while IFS=$'\t' read -r scaffold fragment start stop; do
    # Skip lines that do not have all required fields
    if [ -z "$scaffold" ] || [ -z "$fragment" ] || [ -z "$start" ] || [ -z "$stop" ]; then
        continue
    fi

    echo "Processing $scaffold $fragment from $start to $stop"

    # Calculate mean and median coverage for each prophage region
    prophage_mean=$(calculate_mean_coverage "$scaffold" "$start" "$stop")
    prophage_median=$(calculate_median_coverage "$scaffold" "$start" "$stop")
    echo "Prophage region mean coverage: $prophage_mean"
    echo "Prophage region median coverage: $prophage_median"

    # Get the total length of the scaffold from the BAM header
    scaffold_length=$(samtools view -H "$SORTED_BAM_FILE" | grep "@SQ" | grep "SN:${scaffold}" | cut -f3 | cut -d':' -f2)
    
    if [ -z "$scaffold_length" ]; then
        echo "Error: Could not determine the length of scaffold $scaffold."
        continue
    fi

    # Calculate mean and median coverage for the entire scaffold
    entire_mean=$(calculate_mean_coverage "$scaffold" 1 "$scaffold_length")
    entire_median=$(calculate_median_coverage "$scaffold" 1 "$scaffold_length")
    echo "Entire scaffold mean coverage: $entire_mean"
    echo "Entire scaffold median coverage: $entire_median"

    # Create a temporary file to store non-prophage regions
    tmp_non_prophage=$(mktemp)
    
    # Collect all prophage regions for the scaffold and merge overlapping/adjacent regions
    tail -n +2 "$COORD_FILE" | awk -v scaffold="$scaffold" -F'\t' '$1 == scaffold {print $3, $4}' | sort -k1,1n | \
    awk 'BEGIN {OFS="\t"} {if (NR == 1) {start=$1; end=$2} else {if ($1 <= end+1) {end=($2 > end ? $2 : end)} else {print start, end; start=$1; end=$2}}} END {print start, end}' | \
    awk -v scaffold_length="$scaffold_length" 'BEGIN {OFS="\t"; prev_end=0} {if ($1 > prev_end + 1) print prev_end + 1, $1 - 1; prev_end = $2} END {if (prev_end < scaffold_length) print prev_end + 1, scaffold_length}' > "$tmp_non_prophage"

    # Calculate mean and median coverage for non-prophage regions
    non_prophage_mean=0
    non_prophage_median=0
    non_prophage_count=0

    while IFS=$'\t' read -r np_start np_end; do
        if [ "$np_start" -le "$np_end" ]; then
            np_mean=$(calculate_mean_coverage "$scaffold" "$np_start" "$np_end")
            np_median=$(calculate_median_coverage "$scaffold" "$np_start" "$np_end")
            non_prophage_mean=$(awk "BEGIN {print ($non_prophage_mean * $non_prophage_count + $np_mean * ($np_end - $np_start + 1)) / ($non_prophage_count + ($np_end - $np_start + 1))}")
            non_prophage_median=$(awk "BEGIN {print ($non_prophage_median * $non_prophage_count + $np_median * ($np_end - $np_start + 1)) / ($non_prophage_count + ($np_end - $np_start + 1))}")
            non_prophage_count=$((non_prophage_count + np_end - np_start + 1))
        fi
    done < "$tmp_non_prophage"

    rm -f "$tmp_non_prophage"

    # Calculate the ratio of phage to host coverage for mean and median
    if [ "$non_prophage_mean" != 0 ]; then
        mean_coverage_ratio=$(awk "BEGIN {print $prophage_mean / $non_prophage_mean}")
    else
        mean_coverage_ratio="Inf"
    fi

    if [ "$non_prophage_median" != 0 ]; then
        median_coverage_ratio=$(awk "BEGIN {print $prophage_median / $non_prophage_median}")
    else
        median_coverage_ratio="Inf"
    fi

    echo "Mean coverage ratio (phage/host): $mean_coverage_ratio"
    echo "Median coverage ratio (phage/host): $median_coverage_ratio"

    # Write the results to the output file
    echo -e "${scaffold}\t${fragment}\t${prophage_mean}\t${prophage_median}\t${entire_mean}\t${entire_median}\t${non_prophage_mean}\t${non_prophage_median}\t${mean_coverage_ratio}\t${median_coverage_ratio}" >> "$OUTPUT_FILE"

    echo "-----------------------------"
done
