#!/bin/bash

# ---------------------------------------------------------------
# ready_for_MSA.sh  —  Fully patched safe version (HiTE / panHiTE)
# ---------------------------------------------------------------

if [ $# -ne 3 ]
then
    echo -e "\nUSAGE: $(basename $0) <input.multi.fa> <max> <keep_largest>"
    echo -e "\nCreates a multifasta subset containing <max> sequences:"
    echo -e "  - Takes <keep_largest> longest sequences"
    echo -e "  - Takes (max - keep_largest) random sequences"
    echo -e "\nOutputs:"
    echo -e "  <input>.rdmSubset.fa   (subset fasta)"
    echo -e "  <input>.rdmSubset.len  (length index)"
    exit
fi

fasta="$1"
max=$2
largest=$3

samtools faidx "$fasta"
if [ ! -f "$fasta.fai" ]; then
    echo "[ERROR] Failed to index $fasta"
    exit 1
fi

# Number of sequences in the input FASTA
temp=$(wc -l < "$fasta.fai")

# ---------------------------------------------------------------
# Case 1: number of sequences <= max  → passthrough
# ---------------------------------------------------------------
if [ "$temp" -le "$max" ]
then
    echo "[INFO] $fasta has only $temp sequences (<= $max). No subsampling needed."

    # Clean previous outputs
    : > "$fasta.rdmSubset.fa"
    : > "$fasta.rdmSubset.len"

    # Output entire FASTA
    cp "$fasta" "$fasta.rdmSubset.fa"

    # Create .len file (sorted by length)
    awk '{print $1,$2}' "$fasta.fai" | sort -nk 2 -r > "$fasta.rdmSubset.len"

    rm "$fasta.fai"
    exit 0
fi

# ---------------------------------------------------------------
# Case 2: perform largest + random sampling
# ---------------------------------------------------------------
if [ "$largest" -lt 0 ] || [ "$largest" -gt "$max" ]; then
    echo "[ERROR] <keep_largest> must be between 0 and <max>"
    rm "$fasta.fai"
    exit 1
fi

# Sort fasta.fai by length descending
sort -nk 2 -r "$fasta.fai" > "$fasta.fai.sortedbylength"

# Number to sample randomly = max - largest
r=$((max - largest))

# Clean temporary output files
: > "$fasta.headers"
: > "$fasta.rdmSubset.len"

# Step 1: pick the largest sequences
head -n "$largest" "$fasta.fai.sortedbylength" | awk '{print $1}' >> "$fasta.headers"

# Step 2: pick r random sequences (avoid head -0)
awk "NR > $largest {print}" "$fasta.fai.sortedbylength" > "$fasta.temp"

if [ "$r" -gt 0 ]; then
    sort -R "$fasta.temp" | head -n "$r" | awk '{print $1}' >> "$fasta.headers"
fi

# Step 3: extract the sequences using grep -a (panSN safe)
grep -a -A 1 -f "$fasta.headers" "$fasta" | sed '/^--/d' > "$fasta.rdmSubset.fa"

# Step 4: produce len file sorted by length desc
awk '
    FNR==NR {a[$1]=$2; next}
    {print $1, a[$1]}
' "$fasta.fai" "$fasta.headers" | sort -nk 2 -r > "$fasta.rdmSubset.len"

# Clean temp files
rm "$fasta.fai.sortedbylength" "$fasta.headers" "$fasta.temp" "$fasta.fai"

echo "[INFO] Subset written to:"
echo "       $fasta.rdmSubset.fa"
echo "       $fasta.rdmSubset.len"
