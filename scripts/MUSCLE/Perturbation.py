import random
from Bio.Align import MultipleSeqAlignment


def perturb_alignment(alignment):
    if alignment is None or len(alignment) <= 2:
        return alignment, []

    num_to_remove = max(1, int(0.1 * len(alignment)))

    alignment_list = list(alignment)
    removed = random.sample(alignment_list, num_to_remove)
    remaining = MultipleSeqAlignment([seq for seq in alignment_list if seq not in removed])

    return remaining, removed
