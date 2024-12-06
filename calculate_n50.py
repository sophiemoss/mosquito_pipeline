from Bio import SeqIO
import argparse

def calculate_n50(fasta_file):
    # Parse the FASTA file and get sequence lengths
    sequence_lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

    # Sort the sequence lengths in descending order
    sorted_lengths = sorted(sequence_lengths, reverse=True)

    # Calculate the total length of all sequences
    total_length = sum(sorted_lengths)

    # Calculate the N50
    cumulative_length = 0
    n50 = 0

    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            n50 = length
            break

    return n50

def main():
    parser = argparse.ArgumentParser(description="Calculate N50 from a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")

    args = parser.parse_args()

    try:
        n50_value = calculate_n50(args.fasta_file)
        print(f"N50 value: {n50_value}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

