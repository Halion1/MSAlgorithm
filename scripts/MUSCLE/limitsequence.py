from Bio import SeqIO


def limit_sequences(input_fasta, output_fasta, num_sequences):
    """
    Read sequences from input_fasta, keep only the first num_sequences,
    and write them to output_fasta.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        # Parse the sequences from the input FASTA file
        sequences = list(SeqIO.parse(infile, 'fasta'))

        # Keep only the first num_sequences
        limited_sequences = sequences[:num_sequences]

        # Write the limited set of sequences to the output FASTA file
        SeqIO.write(limited_sequences, outfile, 'fasta')