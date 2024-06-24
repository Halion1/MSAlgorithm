import math

def kimura_distance(seq1, seq2):
    # Ensure the sequences are of equal length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")

    transitions = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
    transversions = [("A", "C"), ("C", "A"), ("G", "T"), ("T", "G"),
                     ("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")]

    ts_count = 0
    tv_count = 0
    for n1, n2 in zip(seq1, seq2):
        if (n1, n2) in transitions:
            ts_count += 1
        elif (n1, n2) in transversions:
            tv_count += 1

    p = ts_count / len(seq1)
    q = tv_count / len(seq1)

    try:
        K = -0.5 * math.log(1 - 2*p - q) - 0.25 * math.log(1 - 2*q)
        return K
    except ValueError:
        # This will occur if the argument of the logarithm is non-positive.
        # It indicates an excessive sequence divergence.
        return float('inf')  # Infinite distance