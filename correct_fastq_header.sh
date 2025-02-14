#!/bin/bash

filename=$1

# Check if filename is provided
if [ -z "$filename" ]; then
    echo "Please provide the path to the FASTQ file as an argument."
    exit 1
fi

# Check if the file exists
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi

# Create a temporary file
tmpfile=$(mktemp)

# Correct the header lines and write to the temporary file
awk '{if ($1 ~ /^@ERR/) {split($1, a, " "); $1 = "@" a[1] " " a[2]} print}' "$filename" > "$tmpfile"

# Replace the original file with the corrected lines
mv "$tmpfile" "$filename"

echo "Header lines corrected successfully."
