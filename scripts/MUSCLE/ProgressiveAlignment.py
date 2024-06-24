from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import MultipleSeqAlignment
from pad import pad_sequences

def progressive_alignment(guide_tree, sequences):
    sequences_dict = {(seq.id,): seq for seq in sequences}
    alignments = dict(sequences_dict)  # Shallow copy of sequences_dict

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1.5
    aligner.extend_gap_score = -0.5

    for seq in sequences:
        assert (seq.id,) in alignments, f"Missing sequence: {seq.id}"

    def align_profiles(node):
        if node.is_terminal():  # If it's a leaf node
            seq = alignments.get((node.name,))
            if not seq:
                print(f"Error: No sequence found for terminal node: {node.name}")
                return None
            return MultipleSeqAlignment([seq])

        # Recursively get alignments for left and right children
        left_alignment = align_profiles(node.clades[0])
        right_alignment = align_profiles(node.clades[1])

        if not left_alignment or not right_alignment:
            if not left_alignment:
                print(f"Error: No alignment found for left child of node: {node.name}")
            if not right_alignment:
                print(f"Error: No alignment found for right child of node: {node.name}")
            return None

        # Convert MSA to single sequence for left and right
        left_seq = "".join([str(record.seq) for record in left_alignment])
        right_seq = "".join([str(record.seq) for record in right_alignment])

        # Pad sequences so they have the same length
        left_seq, right_seq = pad_sequences(left_seq, right_seq)

        # Create SeqRecords
        left_record = SeqRecord(Seq(left_seq))
        right_record = SeqRecord(Seq(right_seq))

        alignment_pair = aligner.align(left_record.seq, right_record.seq)
        if not alignment_pair:
            print("Error: Pairwise alignment failed")
            return None

        aligned_left = alignment_pair[0][0]
        aligned_right = alignment_pair[0][1]

        # Use tuple as key for the combined alignment
        combined_key = (node.clades[0].name, node.clades[1].name)

        # Ensure we are not overwriting any key
        assert combined_key not in alignments, "Attempting to overwrite an existing key in alignments"

        combined_alignment = MultipleSeqAlignment([
            SeqRecord(Seq(aligned_left), id=combined_key[0]),
            SeqRecord(Seq(aligned_right), id=combined_key[1])
        ])

        alignments[combined_key] = combined_alignment
        return combined_alignment

    final_alignment = align_profiles(guide_tree.root)
    if not final_alignment:
        print("Failed to generate alignment")
        return None

    return final_alignment
