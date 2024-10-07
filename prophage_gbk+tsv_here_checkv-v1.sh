#!/bin/bash

# Check if the user provided the required argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_virsorter_output>"
    exit 1
fi

# Warn the user if BASE_DIR ends with a trailing slash
if [[ "$1" == */ ]]; then
    echo "Warning: The provided path should not end with a trailing slash (/)."
    echo "This may affect proper file operations within the script."
    exit 1
fi

echo "Initializing script with provided path: $1"

# Initialize Conda (assuming this is correct and intentional)
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    echo "Initializing Conda..."
    source ~/miniconda3/etc/profile.d/conda.sh
else
    echo "Conda initialization script not found. Please check the Conda installation path."
    exit 1
fi

# Activate Conda environment with bedtools and biopython
echo "Activating Conda environment 'bedtools_env'..."
conda activate bedtools_env
if [ $? -ne 0 ]; then
    echo "Failed to activate Conda environment 'bedtools_env'. Please check if the environment exists."
    exit 1
fi
echo "Conda environment 'bedtools_env' activated successfully."

# Define the base directory based on the provided argument
BASE_DIR="$1"

# Split the pharokka gbk file into its subfiles
echo "Splitting GenBank file..."
python ./processing_scripts_checkv_v1/split_genbank.py "$BASE_DIR/pharokka_output/pharokka.gbk" "$BASE_DIR/pharokka_output/"
if [ $? -ne 0 ]; then
    echo "Failed to execute split_genbank.py script."
    exit 1
fi
echo "GenBank file split successfully."

# Copy pharokka_*_subset.gbk files to the current directory and rename them with $BASE_DIR prefix
echo "Copying and renaming subset GenBank files..."
for gbk_file in "$BASE_DIR/pharokka_output/pharokka_"*.gbk; do
    if [ -f "$gbk_file" ]; then
        filename=$(basename "$gbk_file")
        cp "$gbk_file" "./${BASE_DIR}_$filename"
        echo "Copied and renamed: $filename"
    else
        echo "No subset gbk files found in $BASE_DIR/pharokka_output/"
    fi
done

# Copy tsv files to current directory and rename them with $BASE_DIR prefix
echo "Copying and renaming TSV files..."
if [ -f "$BASE_DIR/script_output.tsv" ]; then
    cp "$BASE_DIR/script_output.tsv" "./${BASE_DIR}_script_output.tsv"
    echo "Copied and renamed: script_output.tsv"
    cp "$BASE_DIR/manual_prophage_coordinates.tsv" "./${BASE_DIR}_manual_prophage_coordinates.tsv"
    echo "Copied and renamed: manual_prophage_coordinates.tsv"
else
    echo "script_output.tsv not found in $BASE_DIR."
    exit 1
fi

echo "Script execution completed successfully."
