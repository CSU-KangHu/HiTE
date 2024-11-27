#!/usr/bin/env python3
import os
import re
import subprocess
import argparse
import sys
import concurrent.futures
import tempfile
import shutil
import pandas as pd

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cur_dir)


def run_command(cmd):
    """
    Run a shell command and check for errors.
    """
    print(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        raise RuntimeError(result.stderr)
    return result.stdout


def convert_gff_to_gtf(input_gff, output_gtf, temp_dir):
    """
    Convert a GFF file to Ensembl-compatible GTF format using gffread and gtftk.
    """
    print("Converting GFF to GTF using gffread...")
    temp_gtf = os.path.join(temp_dir, "temp_converted.gtf")
    run_command(f"gffread {input_gff} -T -o {temp_gtf}")

    print("Ensuring GTF is in Ensembl-compatible format using gtftk...")
    run_command(f"gtftk convert_ensembl -i {temp_gtf} > {output_gtf}")

    # Remove intermediate file
    os.remove(temp_gtf)
    print(f"Conversion complete. Output GTF: {output_gtf}")


# def extract_coding_genes(input_gtf, temp_dir):
#     """
#     Extract coding genes directly to a temporary GTF file.
#     """
#     temp_coding_mapid = os.path.join(temp_dir, "temp_coding_mapid.txt")
#     temp_coding_gtf = os.path.join(temp_dir, "temp_coding.gtf")
#     temp_cds_gtf = os.path.join(temp_dir, "temp_cds.gtf")
#
#     run_command(f"gtftk select_by_key -k feature -v CDS -i {input_gtf} > {temp_cds_gtf}")
#     run_command(f"gtftk tabulate -k transcript_id,gene_id --unique --no-header -i {temp_cds_gtf} > {temp_coding_mapid}")
#     run_command(f"perl " + cur_dir + f"/gtf_extract.pl {input_gtf} {temp_coding_mapid} > {temp_coding_gtf}")
#
#     return temp_coding_gtf


def extract_longest_transcripts(input_gtf, output_gtf):
    """
    Extract longest transcripts from the input GTF file to the output GTF.
    """
    run_command(f"gtftk short_long -g -l -i {input_gtf} > {output_gtf}")


# Parse attributes to extract transcript ID and gene ID
def parse_attributes(attr_str):
    print(attr_str)
    attrs = {}
    for item in attr_str.strip(';').split('; '):
        if ' ' in item:
            key, value = item.split(' ', 1)
            attrs[key] = value.strip('"')  # 移除引号
    return attrs

def extract_longest_transcripts(input_gtf, output_gtf):
    """
    Extract longest transcripts from the input GTF file to the output GTF.
    This implementation reads the GTF file, identifies the longest transcript for each gene,
    and writes it back to a new GTF file.
    """
    # Read GTF file into a pandas DataFrame
    column_names = [
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ]
    gtf = pd.read_csv(input_gtf, sep='\t', comment='#', names=column_names, low_memory=False)

    # Filter for transcript features
    transcripts = gtf[gtf["feature"] == "transcript"]

    transcripts[["transcript_id", "gene_id"]] = transcripts["attribute"].apply(parse_attributes).apply(pd.Series)

    # Calculate transcript lengths and find the longest for each gene
    transcripts["length"] = transcripts["end"] - transcripts["start"]
    longest_transcripts = transcripts.loc[transcripts.groupby("gene_id")["length"].idxmax()]

    # Filter the original GTF file to include only the longest transcripts
    longest_ids = set(longest_transcripts["transcript_id"])
    longest_gtf = gtf[gtf["attribute"].str.contains('|'.join(longest_ids))]

    # Write to the output GTF file
    longest_gtf.to_csv(output_gtf, sep='\t', header=False, index=False)


def process_gtf(input_file, output_dir):
    """
    Process a single GTF/GFF file: Convert GFF to GTF if needed, extract coding genes, and extract longest transcripts.
    """
    base_name = os.path.basename(input_file)
    temp_output_gtf = os.path.join(output_dir, os.path.splitext(base_name)[0] + ".gtf")

    # Create a unique temporary directory for each file processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Determine file type and convert GFF to GTF if necessary
        if input_file.endswith(".gff") or input_file.endswith(".gff3"):
            temp_gtf = os.path.join(temp_dir, os.path.splitext(base_name)[0] + ".converted.gtf")
            convert_gff_to_gtf(input_file, temp_gtf, temp_dir)
            input_gtf = temp_gtf
        elif input_file.endswith(".gtf"):
            input_gtf = input_file
        else:
            raise ValueError("Unsupported file format. Please provide a GFF or GTF file.")

        # Step 1: Extract coding genes to a temporary GTF
        coding_gtf = extract_coding_genes(input_gtf, temp_dir)

        # Step 2: Extract longest transcripts and save to the output directory
        extract_longest_transcripts(coding_gtf, temp_output_gtf)

        print(f"Processing completed for {input_file}. Output GTF: {temp_output_gtf}")



def main():
    parser = argparse.ArgumentParser(
        description="Process multiple GFF/GTF files to extract coding genes and longest transcripts.")
    parser.add_argument("--input_files", required=True, nargs='+', help="List of input GFF/GTF files.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Use ThreadPoolExecutor to run the processing in parallel for each input file
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Map input files to the processing function
        futures = [executor.submit(process_gtf, input_file, args.output_dir) for input_file in args.input_files]

        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            future.result()  # Raise any exceptions that occurred during the execution

    print("All files processed successfully.")


if __name__ == "__main__":
    main()
