#!/usr/bin/env python
import gzip
from Bio import SeqIO
import os

def get_barcode_from_header(header):
    parts = header.split()
    for part in parts:
        if part.startswith("barcode="):
            return part.split('=')[1]
    return None

def split_fastq_by_barcode(fastq_file):
    barcodes = {}

  # opens g zip file
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            barcode = get_barcode_from_header(record.description)
            
            if barcode:
                if barcode not in barcodes:
                    barcodes[barcode] = []
                barcodes[barcode].append(record)

    for barcode, records in barcodes.items():
        output_filename = f"{barcode}.fastq"
        with open(output_filename, "w") as output_handle:
            SeqIO.write(records, output_handle, "fastq")
# Usage example
split_fastq_by_barcode("file.fastq.gz") 
