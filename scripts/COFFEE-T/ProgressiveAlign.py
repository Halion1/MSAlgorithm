from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from PairwiseAlignmentScore import needleman_wunsch



def align_sequences(seq1, seq2, substitution_matrix, gap_penalty):
    # Call Needleman-Wunsch to get the score and the aligned sequences.
    _, align1, align2 = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
    return align1, align2

def align_sequence_groups(group1, group2, substitution_matrix, gap_penalty):
    # Calculate the consensus sequence for each group
    consensus_group1 = calculate_consensus(group1)
    consensus_group2 = calculate_consensus(group2)

    # Align the consensus sequences
    align1, align2 = align_sequences(consensus_group1, consensus_group2, substitution_matrix, gap_penalty)

    # Insert gaps into the original sequences based on the consensus alignment
    aligned_group1 = [insert_gaps(seq, align1) for seq in group1]  # You need to implement insert_gaps
    aligned_group2 = [insert_gaps(seq, align2) for seq in group2]  # Similarly, implement insert_gaps

    # Combine the groups into a new aligned group
    new_group = aligned_group1 + aligned_group2

    return new_group


def insert_gaps(seq, consensus):
    print("Original sequence:", seq)
    print("Consensus:", consensus)
    aligned_seq = ""
    seq_index = 0

    for i in range(len(consensus)):
        if consensus[i] == '-':  # if there's a gap in the consensus sequence
            aligned_seq += '-'  # insert a gap in the original sequence
        else:  # if there's no gap
            if seq_index < len(seq):
                aligned_seq += seq[seq_index]  # copy the character from the original sequence
                seq_index += 1  # move to the next character in the original sequence
            else:
                aligned_seq += '-'  # if the original sequence is shorter, insert a gap

    # if the original sequence is longer, append the remaining parts to the aligned sequence
    if seq_index < len(seq):
        aligned_seq += seq[seq_index:]

    return aligned_seq


def traverse_tree(node, substitution_matrix, gap_penalty, sequence_dict):
    # If the sequence for the current node is already computed, return it
    if node.name in sequence_dict:
        return sequence_dict[node.name]

    # If the node is a leaf, return its sequence
    if node.is_terminal():
        return sequence_dict[node.name]

    # Recursive alignment of child sequences
    aligned_sequences = []
    for child in node.clades:
        aligned_seq = traverse_tree(child, substitution_matrix, gap_penalty, sequence_dict)
        aligned_sequences.append(aligned_seq)

    # Align each child sequence with the previous consensus
    consensus_sequence = aligned_sequences[0]  # start with the first sequence
    for aligned_seq in aligned_sequences[1:]:
        consensus_sequence, aligned_seq = align_sequences(consensus_sequence, aligned_seq, substitution_matrix, gap_penalty)

    # Update the sequence dictionary with the computed consensus
    sequence_dict[node.name] = consensus_sequence

    return consensus_sequence

def create_consensus(seq1, seq2):
    consensus = []
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            consensus.append(s1)
        elif s1 == '-' or s2 == '-':
            consensus.append(s1 if s1 != '-' else s2)
        else:
            consensus.append(s1)  # selecting from the first sequence
    return ''.join(consensus)


def progressive_alignment(guide_tree, substitution_matrix, gap_penalty, sequences, names=None):
    """
    Perform progressive alignment of sequences using a guide tree.

    :param guide_tree: The guide tree (Bio.Phylo.BaseTree.Tree) used for alignment.
    :param substitution_matrix: Substitution matrix used for scoring alignments.
    :param gap_penalty: Penalty for introducing a gap in the alignment.
    :param sequences: List of sequences (strings) to be aligned.
    :param names: List of original names of sequences.
    :return: MultipleSeqAlignment object containing the final alignment.
    """

    # If original names aren't provided, default to Seq1, Seq2, etc.
    if names is None:
        names = [f"Seq{i + 1}" for i in range(len(sequences))]

    # Map sequence names to actual sequences
    sequence_dict = {name: sequences[i] for i, name in enumerate(names)}
    print("Guide tree:", guide_tree)
    print("Original sequences:", sequences)

    # Traverse the tree and get the aligned sequences
    aligned_sequences = traverse_tree(guide_tree.root, substitution_matrix, gap_penalty, sequence_dict)

    # Convert the list of aligned sequences to a MultipleSeqAlignment object
    msa = MultipleSeqAlignment(
        [SeqRecord(Seq(seq), id=name) for name, seq in zip(names, aligned_sequences)]
    )
    print("Guide tree:", guide_tree)
    print("Original sequences:", sequences)

    return msa


def calculate_consensus(sequences):
    consensus_sequence = ""
    sequence_length = len(sequences[0])  # Assumes all sequences are aligned and have the same length

    for i in range(sequence_length):
        # Count occurrences of each character at position i
        char_count = {}
        for seq in sequences:
            char = seq[i]
            if char in char_count:
                char_count[char] += 1
            else:
                char_count[char] = 1

        # Find the character with the maximum count
        consensus_char = max(char_count, key=char_count.get)
        consensus_sequence += consensus_char

    return consensus_sequence


def check_sequences(sequences):
    # Function that checks if all sequences are the same length
    length = len(sequences[0])
    for seq in sequences:
        if len(seq) != length:
            raise ValueError("Sequences are not all the same length!")
