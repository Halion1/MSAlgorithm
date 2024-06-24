def get_consensus_dna(aligned_sequences):
    consensus = []

    # Iterate over each column in the alignment
    for col in zip(*aligned_sequences):
        # Count each nucleotide in the column
        nucleotide_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for nucleotide in col:
            if nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] += 1

        # Find the nucleotide with the highest count for this column
        most_common_nucleotide = max(nucleotide_counts, key=nucleotide_counts.get)

        # If two nucleotides have the same count, we append 'N'
        max_count = nucleotide_counts[most_common_nucleotide]
        if list(nucleotide_counts.values()).count(max_count) > 1:
            consensus.append('N')
        else:
            consensus.append(most_common_nucleotide)

    return ''.join(consensus)
