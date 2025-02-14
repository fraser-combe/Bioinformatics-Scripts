import sys

sam_file = sys.argv[1]
reference_file = sys.argv[2]
output_file = sys.argv[3]

# Parse reference lengths
ref_lengths = {}
with open(reference_file, 'r') as ref:
    ref_name = None
    for line in ref:
        if line.startswith('>'):
            ref_name = line[1:].strip()
            ref_lengths[ref_name] = 0
        else:
            ref_lengths[ref_name] += len(line.strip())

# Filter SAM file
with open(sam_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('@'):
            outfile.write(line)
            continue
        fields = line.strip().split('\t')
        flag = int(fields[1])
        ref_name = fields[2]
        read_length = len(fields[9])

        if flag & 4 == 0:
            ref_length = ref_lengths.get(ref_name, 0)
            if read_length >= 0.99 * ref_length:
                outfile.write(line)
