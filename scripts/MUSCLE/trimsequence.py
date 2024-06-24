from Bio import SeqIO


def trim_sequences(input_fasta, output_fasta, desired_length):
    """
    Read sequences from input_fasta, trim them to desired_length,
    and write the trimmed sequences to output_fasta.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        # Parse the sequences from the input FASTA file
        sequences = list(SeqIO.parse(infile, 'fasta'))

        # Trim each sequence to the desired length
        for seq in sequences:
            seq.seq = seq.seq[:desired_length]

        # Write the trimmed sequences to the output FASTA file
        SeqIO.write(sequences, outfile, 'fasta')


