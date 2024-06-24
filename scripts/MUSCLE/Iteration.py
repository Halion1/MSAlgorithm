from AlignmentScoring import sum_of_pairs_score
from Perturbation import perturb_alignment
from ProgressiveAlignment import progressive_alignment

def iterative_refinement(initial_alignment, guide_tree, num_iterations=10):
    current_alignment = initial_alignment
    current_score = sum_of_pairs_score(current_alignment)

    print(f"Initial alignment contains: {[record.id for record in initial_alignment]}")

    for i in range(num_iterations):
        remaining, removed = perturb_alignment(current_alignment)

        print(f"Remaining after perturbation: {[seq.id for seq in remaining]}")
        print(f"Removed during perturbation: {[seq.id for seq in removed]}")

        combined_sequences = list(remaining) + list(removed)
        new_alignment = progressive_alignment(guide_tree, combined_sequences)
        new_score = sum_of_pairs_score(new_alignment)

        if new_score > current_score:
            current_alignment = new_alignment
            current_score = new_score

    return current_alignment


