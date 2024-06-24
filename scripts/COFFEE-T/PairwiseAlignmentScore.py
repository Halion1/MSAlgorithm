import numpy as np
def create_substitution_matrix(match_score, mismatch_score, characters):
    """
    Create a substitution matrix given a set of characters and scores for match and mismatch.
    Gap penalties are handled separately.
    """
    substitution_matrix = {}
    for char1 in characters:
        for char2 in characters:
            if char1 == char2:
                substitution_matrix[(char1, char2)] = match_score
            else:
                substitution_matrix[(char1, char2)] = mismatch_score

    # Add gap penalties
    gap_score = mismatch_score  # or any other value you deem appropriate
    for char in characters:
        substitution_matrix[('-', char)] = gap_score
        substitution_matrix[(char, '-')] = gap_score
    substitution_matrix[('-', '-')] = gap_score  # typically the same as the gap extension penalty

    return substitution_matrix


def needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty):
    """
    Implement the Needleman-Wunsch algorithm for global sequence alignment.
    """
    # Initialize the score matrix.
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)

    # Initialization: penalize gaps at the beginning of sequences.
    for i in range(1, rows):
        score_matrix[i][0] = i * gap_penalty
    for j in range(1, cols):
        score_matrix[0][j] = j * gap_penalty

    # Fill the score matrix and traceback matrix.
    traceback_matrix = np.zeros((rows, cols), dtype=int)  # 0: diagonal, 1: up, 2: left
    for i in range(1, rows):
        for j in range(1, cols):
            # Lookup score in the substitution matrix; default to gap_penalty if not found.
            match = score_matrix[i - 1][j - 1] + substitution_matrix.get((seq1[i - 1], seq2[j - 1]), gap_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            # Find the action with the highest score.
            actions = [match, delete, insert]
            best_score = max(actions)
            best_action = actions.index(best_score)

            score_matrix[i][j] = best_score
            traceback_matrix[i][j] = best_action

    # Traceback: reconstruct the alignment from the matrices.
    align1, align2 = [], []
    i, j = rows - 1, cols - 1
    while i > 0 or j > 0:
        action = traceback_matrix[i][j]
        if action == 0:  # Diagonal move
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif action == 1:  # Move up
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1
        else:  # Move left
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1

    # Reverse alignments to get the correct order.
    align1 = ''.join(align1[::-1])
    align2 = ''.join(align2[::-1])

    return score_matrix[rows - 1][cols - 1], align1, align2

# Keeping the align_sequences function unchanged
def align_sequences(seq1, seq2, substitution_matrix, gap_penalty):
    _, align1, align2 = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
    return align1, align2


def pairwise_alignment_score(seq1, seq2, characters):
    """
    Calculate the pairwise alignment score between two sequences.
    """
    # Define the match and mismatch scores.
    match_score = 2
    mismatch_score = -1
    gap_penalty = -1

    # Include gap character as a valid character for the substitution matrix.
    characters = set(characters)
    characters.add('-')

    # Create a substitution matrix.
    substitution_matrix = create_substitution_matrix(match_score, mismatch_score, characters)

    # There's no need to sanitize sequences by removing gaps as they're now valid characters.
    # clean_seq1 = sanitize_sequence(seq1)
    # clean_seq2 = sanitize_sequence(seq2)

    # Compute the alignment score.
    score = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
    return score


def sanitize_sequence(seq):
    """
    Remove gap characters from the sequence.
    """
    return seq.replace('-', '')
