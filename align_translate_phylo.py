from collections import Counter
import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

def write_consensus_to_fasta(consensus, output_fasta, output_protein_fasta):
    nucleotide_records = []
    protein_records = []
    for name, seq in consensus.items():
        # Create nucleotide SeqRecord
        nucleotide_record = SeqRecord(Seq(seq), id=f"{name}_nucleotide")
        nucleotide_records.append(nucleotide_record)
        
        # Translate to amino acid sequence
        protein_seq = Seq(seq).translate()
        
        # Create protein SeqRecord
        protein_record = SeqRecord(protein_seq, id=f"{name}_protein")
        protein_records.append(protein_record)

    # Write nucleotide sequences to fasta
    SeqIO.write(nucleotide_records, output_fasta, "fasta")
    
    # Write protein sequences to fasta
    SeqIO.write(protein_records, output_protein_fasta, "fasta")
def create_consensus(aligned_sequences, threshold=0.0):
    consensus = {}
    for name, seqs in aligned_sequences.items():
        consensus_seq = []
        min_length = min(len(seq) for seq in seqs if len(seq) > 0)
        
        if min_length == 0:
            print(f"Warning: Skipped {name} due to zero-length sequence")
            continue

        for i in range(min_length):
            bases_at_i = [seq[i] for seq in seqs if len(seq) > i]
            total_count = len(bases_at_i)
            base_count = Counter(bases_at_i)
            
            # Calculate fractions and find those above the threshold
            above_threshold = {base: count / total_count for base, count in base_count.items() if (count / total_count) >= threshold}
            
            if above_threshold:
                most_common_base, _ = max(above_threshold.items(), key=lambda x: x[1])
                consensus_seq.append(most_common_base)
            else:
                # Handle ambiguous base, or simply add 'N' if no base meets the threshold
                consensus_seq.append('N')

        consensus[name] = "".join(consensus_seq)
        
    return consensus

# Create a combined initial FASTA with both reference and consensus sequences
def create_combined_fasta(reference_fasta, consensus_fasta, output_fasta):
    # Read the consensus sequences from the provided fasta file
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))

    # Read the reference sequences from the provided fasta file
    reference_records = list(SeqIO.parse(reference_fasta, "fasta"))

    # Combine both lists
    combined_records = reference_records + consensus_records

    # Write to output fasta file
    SeqIO.write(combined_records, output_fasta, "fasta")

def parse_sam(sam_file):
    aligned_sequences = {}
    with open(sam_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            rname = fields[2]  # RNAME field
            seq = fields[9]  # SEQ field

            #print(f"Debug: rname={rname}, seq={seq}")  # Debugging line
            
            if seq == '*':  # Skip lines without sequence data
                continue

            if rname in aligned_sequences:
                aligned_sequences[rname].append(seq)
            else:
                aligned_sequences[rname] = [seq]
    return aligned_sequences

def create_fasta(aligned_sequences, output_fasta):
    records = []
    for name, seqs in aligned_sequences.items():
        combined_seq = "".join(seqs)  # Assuming that you want to concatenate sequences
        records.append(SeqRecord(Seq(combined_seq), id=name))
    SeqIO.write(records, output_fasta, "fasta")

def run_mafft(input_fasta, output_fasta):
    print("Running MAFFT alignment...")
    result = subprocess.run(["mafft", input_fasta], stdout=open(output_fasta, "w"))
    if result.returncode != 0:
        print("Error running MAFFT.")
        exit(1)

#def create_tree(aligned_fasta, output_tree):
 #   root_name = "Root sequence Name"
  #  print("Creating phylogenetic tree...")
   # subprocess.run(["FastTree", "-outgroup", root_name, aligned_fasta], stdout=open(output_tree, "w"))
    
def create_tree_with_iqtree(aligned_fasta, output_tree, outgroup):
    print("Creating phylogenetic tree with IQ-TREE...")
    subprocess.run(["iqtree2", "-s", aligned_fasta, "-o", outgroup, "-nt", "AUTO", "-pre", output_treeIQ])


# Check if files exist
if not os.path.exists("reference.fasta") or not os.path.exists("input.fastq.gz"):
    print("Reference or input files not found.")
    exit(1)

# Align sequences with minimap2
print("Running Minimap2 alignment...")
subprocess.run(["minimap2", "-ax", "map-ont", "-o", "output.sam", "reference.fasta", "input.fastq.gz"])

# Parse SAM file
aligned_sequences = parse_sam("output.sam")

# Check if we have aligned sequences
if not aligned_sequences:
    print("No aligned sequences found.")
    exit(1)

# Create consensus sequences
consensus = create_consensus(aligned_sequences)


# Write the consensus sequences and their translations to fasta files
write_consensus_to_fasta(consensus, "consensus_nucleotide.fasta", "consensus_protein.fasta")


# Create initial FASTA with reference and consensus sequences
initial_fasta = "initial_combined.fasta"
create_combined_fasta("reference.fasta", consensus_fasta, initial_fasta)

# Run MAFFT for alignment
aligned_fasta = "aligned.fasta"
run_mafft(initial_fasta, aligned_fasta)

# Create tree
#output_tree = "tree.newick"
#create_tree(aligned_fasta, output_tree)

# Create tree with IQ-TREE and specify outgroup
outgroup = "Outgroup name"
output_treeIQ = "tree2"
create_tree_with_iqtree(aligned_fasta, output_treeIQ, outgroup)

print("All steps completed successfully!")


