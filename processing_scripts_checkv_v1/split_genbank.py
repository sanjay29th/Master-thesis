import os
import argparse
from Bio import SeqIO

def split_genbank(input_file, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse the input GenBank file
    records = list(SeqIO.parse(input_file, "genbank"))

    # Loop through each record and write to a new file
    for i, record in enumerate(records, 1):
        output_file = os.path.join(output_dir, f"pharokka_{i}.gbk")
        SeqIO.write(record, output_file, "genbank")
        print(f"Written {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a multi-sequence GenBank file into individual files.")
    parser.add_argument("input_file", type=str, help="Path to the input multi-sequence GenBank file.")
    parser.add_argument("output_dir", type=str, help="Directory where the output files will be saved.")
    
    args = parser.parse_args()
    
    split_genbank(args.input_file, args.output_dir)
