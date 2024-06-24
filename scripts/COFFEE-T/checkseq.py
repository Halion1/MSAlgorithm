from Bio.Seq import Seq


def check_sequences(sequences):
    # Function that checks if all sequences are the same length
    length = len(sequences[0])
    for seq in sequences:
        if len(seq) != length:
            raise ValueError(f"Sequences are not all the same length! Found a sequence of length {len(seq)} but expected {length}.")


def pad_sequences(alignment):
    max_length = max(len(record.seq) for record in alignment)
    for record in alignment:
        record.seq = record.seq + Seq('-' * (max_length - len(record.seq)))
    return alignment

