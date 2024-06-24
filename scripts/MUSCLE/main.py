# Import .py files
from kmercounts import compute_distance_matrix
from tree import generate_guide_tree, convert_to_lower_triangle
from ProgressiveAlignment import progressive_alignment
from trimsequence import trim_sequences
from limitsequence import limit_sequences
import time
from Bio import SeqIO

# Example usage:
#trim_sequences("try1.fasta", "try2.fasta", 5000)

# Example usage:
#limit_sequences("try2.fasta", "try3.fasta", 7)



# Read in the sequences
sequences = list(SeqIO.parse("C:/Users/rodri/Desktop/MUSCLE/try3.fasta", "fasta"))

for iteration in range(10):  # Start the loop
    print(f"Iteration {iteration + 1}:")
    start_time = time.time()

    # 2. Compute the initial distance matrix
    matrix = compute_distance_matrix(sequences)
    print(matrix.matrix)  # This will show the distance matrix

    # 3 & 4. Convert and generate guide tree
    sequence_names = [seq.id for seq in sequences]
    assert len(sequence_names) == len(matrix.matrix), "Names and matrix dimensions do not match."

    print("Expected lower triangle matrix:")
    for i in range(1, len(sequence_names)):
        print(matrix.matrix[i][:i])

    biopy_matrix = convert_to_lower_triangle(matrix, sequence_names)
    print("Constructed lower triangle matrix:")
    for row in biopy_matrix:
        print(row)

    guide_tree = generate_guide_tree(matrix, sequence_names)

    # If you want to view the tree uncomment the following line:
    # Phylo.draw_ascii(guide_tree)

    # 5. Generate initial alignment
    initial_alignment = progressive_alignment(guide_tree, sequences)
    if initial_alignment is None:
        print("Failed to generate alignment")
        # Removed the breakpoint() for the loop to continue to the next iteration
        continue

    end_time = time.time()  # Record end time

    print(initial_alignment)
    print(f"Alignment took {end_time - start_time:.2f} seconds.")
    print("------")  # Separator for clarity in output

# Save the alignment to a file, visualize, etc. if needed
# for record in initial_alignment:
#     print(record.id, record.seq)