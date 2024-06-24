def get_consensus_protein(aligned_sequences):
    # Define groups of amino acids based on their properties
    hydrophobic = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'P'}
    polar = {'C', 'G', 'Y', 'N', 'Q', 'S', 'T'} #uncharged
    acidic = {'D', 'E'}
    basic = {'R', 'H', 'K'}

    consensus = []

    for col in zip(*aligned_sequences):
        unique_amino_acids = set(col)

        # If only one amino acid is present, it is the consensus
        if len(unique_amino_acids) == 1:
            consensus.append(unique_amino_acids.pop())
        elif unique_amino_acids.issubset(hydrophobic):
            consensus.append(':')
        elif unique_amino_acids.issubset(polar):
            consensus.append(':')
        elif unique_amino_acids.issubset(acidic):
            consensus.append(':')
        elif unique_amino_acids.issubset(basic):
            consensus.append(':')
        # Here, you can add more cases for other groups as needed
        else:
            consensus.append('.')

    return ''.join(consensus)
