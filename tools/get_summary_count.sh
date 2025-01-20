#!/bin/bash
# Script Function: Process TE annotation files and generate statistical results
# Author: Kang Hu
# Date: 2025-01-19

# Check the number of arguments
if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    echo "Usage: $0 <TE_lib> <HiTE_gff> <genome> [-gff3]"
    echo "  -gff3: Skip steps 1 and 2 if the input is already in GFF3 format."
    exit 1
fi

# Input parameters
TE_lib=$1       # TE library file (FASTA format)
HiTE_gff=$2     # HiTE annotation file (GFF format)
genome=$3       # Genome file (FASTA format)
is_gff3=false   # Flag to check if the input is already in GFF3 format

# Check if the -gff3 flag is provided
if [ $# -eq 4 ] && [ "$4" == "-gff3" ]; then
    is_gff3=true
fi

# Tool paths
gff2RMout="gff2RMout.pl"
count_base="count_base.pl"
buildSummary="buildSummary.pl"

# Step 1: Extract TE family and type information (skip if input is GFF3)
if [ "$is_gff3" = false ]; then
    echo "Extracting TE family and type information..."
    cat "${TE_lib}" | grep '>' | sed 's/>//g' | tr '#' '\t' > panHiTEfam_type.txt
    if [ $? -ne 0 ]; then
        echo "Error: Failed to extract TE family and type information!"
        exit 1
    fi

    # Step 2: Convert GFF file to GFF3 format
    echo "Converting GFF file to GFF3 format..."
    gff2gff3.py panHiTEfam_type.txt "${HiTE_gff}" "${HiTE_gff}.gff3"
    if [ $? -ne 0 ]; then
        echo "Error: GFF file conversion failed!"
        exit 1
    fi
else
    echo "Input is already in GFF3 format. Skipping steps 1 and 2."
    # Use the input GFF file directly as GFF3
    cp "${HiTE_gff}" "${HiTE_gff}.gff3"
fi

# Step 3: Convert GFF3 file to RepeatMasker .out format
echo "Converting GFF3 file to RepeatMasker .out format..."
"${gff2RMout}" "${HiTE_gff}.gff3" "${HiTE_gff}.out"
if [ $? -ne 0 ]; then
    echo "Error: GFF3 file conversion failed!"
    exit 1
fi

# Step 4: Calculate basic genome statistics
echo "Calculating basic genome statistics..."
"${count_base}" "${genome}" > "${genome}.stats"
if [ $? -ne 0 ]; then
    echo "Error: Genome statistics calculation failed!"
    exit 1
fi

# Step 5: Generate TE annotation statistics
echo "Generating TE annotation statistics..."
"${buildSummary}" -maxDiv 40 -stats "${genome}.stats" "${HiTE_gff}.out" > "${HiTE_gff}.tbl" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate statistics!"
    exit 1
fi

# Step 6: Clean up intermediate files
echo "Cleaning up intermediate files..."
if [ "$is_gff3" = false ]; then
    rm -f panHiTEfam_type.txt
fi
rm -f "${HiTE_gff}.gff3" "${HiTE_gff}.out" "${genome}.stats"
if [ $? -ne 0 ]; then
    echo "Warning: Failed to delete some intermediate files!"
fi

echo "Script execution completed! Results saved to: ${HiTE_gff}.tbl"