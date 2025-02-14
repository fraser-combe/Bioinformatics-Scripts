import sys
import pysam

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

# Input files
consensus_fasta = sys.argv[1]
sorted_bam = sys.argv[2]
output_file = sys.argv[3]

# Create a dictionary to hold the read counts for each reference
read_count_dict = {}

# Open the BAM file and count the reads mapped to each reference
with pysam.AlignmentFile(sorted_bam, "rb") as bamfile:
    for read in bamfile:
        if not read.is_unmapped:
            ref_name = bamfile.get_reference_name(read.reference_id)
            read_count_dict[ref_name] = read_count_dict.get(ref_name, 0) + 1

# Read the consensus FASTA file and create new headers
new_lines = []
with open_file(consensus_fasta) as f:
    for line in f:
        if line.startswith('>'):
            ref_name = line[1:].strip()
            read_count = read_count_dict.get(ref_name, 0)
            new_header = f">Sequence_assembled_to_{ref_name}_using_minimap2_{read_count}_reads\n"
            new_lines.append(new_header)
        else:
            new_lines.append(line)

# Write the new FASTA file with renamed headers
with open(output_file, 'w') as f:
    f.writelines(new_lines)
