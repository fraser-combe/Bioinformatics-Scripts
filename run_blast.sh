#!/bin/bash

fasta_file=$1
output_file=$2
error_file=$3

# Run the Perl script with the FASTA file
perl /app/scripts/web_blast.pl $fasta_file > $output_file 2> $error_file
