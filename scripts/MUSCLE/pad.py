def pad_sequences(seq1, seq2):
    """Pad sequences with gap characters to be of the same length."""
    len_diff = len(seq1) - len(seq2)

    if len_diff > 0:
        seq2 += '-' * len_diff
    elif len_diff < 0:
        seq1 += '-' * abs(len_diff)

    return seq1, seq2