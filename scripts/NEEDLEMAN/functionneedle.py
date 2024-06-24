def needlemanwunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)
    matrix = [[0 for j in range(m + 1)] for i in range(n + 1)]

    for i in range(n + 1):
        matrix[i][0] = i * gap
    for j in range(m + 1):
        matrix[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = matrix[i - 1][j] + gap
            insert = matrix[i][j - 1] + gap
            matrix[i][j] = max(match_score, delete, insert)

    align1, align2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (
        match if seq1[i - 1] == seq2[j - 1] else mismatch):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        else:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    return align1[::-1], align2[::-1]


def most_frequent_alignment(alignments):
    """
    Select the most frequent alignment for each sequence.
    """
    from collections import Counter

    # Count occurrences of each alignment
    alignment_counts = Counter(alignments)

    # Return the most common alignment
    return alignment_counts.most_common(1)[0][0]


def multiple_sequence_alignment(sequences):
    if len(sequences) < 2:
        return sequences

    align1, align2 = needlemanwunsch(sequences[0], sequences[1])

    # Dictionary to hold aligned sequences
    aligned_sequences_dict = {sequences[0]: align1, sequences[1]: align2}

    for seq in sequences[2:]:
        new_aligned_sequences_dict = {}
        for previous_seq, aligned_seq in aligned_sequences_dict.items():
            a, b = needlemanwunsch(aligned_seq, seq)
            new_aligned_sequences_dict[previous_seq] = a
            new_aligned_sequences_dict[seq] = b
        aligned_sequences_dict = new_aligned_sequences_dict

    # Create a list of aligned sequences based on the original order
    aligned_sequences = [aligned_sequences_dict[seq] for seq in sequences]

    return aligned_sequences
