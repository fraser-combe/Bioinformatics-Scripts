#!/bin/bash

# Directory containing the FASTA files
FASTA_DIR="path/to/files"

# Find all FASTA files and index them
find "$FASTA_DIR" -type f -name "*.fasta" -or -name "*.fa" | while read fasta_file; do
    echo "Indexing $fasta_file"
    samtools faidx "$fasta_file"
done
