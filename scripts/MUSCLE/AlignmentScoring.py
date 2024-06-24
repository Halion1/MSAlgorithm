from Bio.Align import substitution_matrices


def sum_of_pairs_score(alignment):
    matrix = substitution_matrices.load("BLOSUM62")
    score = 0
    for i in range(len(alignment)):
        for j in range(i+1, len(alignment)):
            for k in range(len(alignment[0])):
                pair = (alignment[i][k], alignment[j][k])
                if '-' in pair:  # gap penalty
                    score -= 1
                else:
                    score += matrix.get(tuple(sorted(pair)), 0)  # use the matrix score
    return score
