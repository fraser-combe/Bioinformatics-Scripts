import matplotlib.pyplot as plt
import sys


def plot_read_length_distribution(fastq_file, sample_name):
    read_lengths = []
    with open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence lines in FASTQ
                read_lengths.append(len(line.strip()))

    if not read_lengths:
        print("No reads found in the input FASTQ file.")
        return

    output_plot = f"{sample_name}_read_length_distribution.png"

    plt.figure(figsize=(10, 6))
    plt.hist(read_lengths, bins=range(min(read_lengths), max(read_lengths) + 1, 10), edgecolor='black')
    plt.title(f"Read Length Distribution for {sample_name}")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Frequency")
    plt.savefig(output_plot)
    plt.close()
    print(f"Plot saved as {output_plot}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python plot_read_length_distribution.py <filtered_fastq> <sample_name>")
        sys.exit(1)

    fastq_file = sys.argv[1]
    sample_name = sys.argv[2]

    plot_read_length_distribution(fastq_file, sample_name)
