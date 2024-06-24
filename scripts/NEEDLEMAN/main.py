import time
from functionneedle import *
from protein import *
from dna import *
from memory_profiler import memory_usage

def read_fasta(file_path):
    """Read sequences from a FASTA file and return as a list."""
    # C:\Users\USUARIO\PycharmProjects\Alingment\spike.txt
    try:
        with open(file_path, 'r') as file:
            sequences = []
            sequence = ''
            for line in file:
                if line.startswith('>'):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ''
                    continue
                sequence += line.strip()
            if sequence:  # To handle the last sequence in the file
                sequences.append(sequence)
        return sequences
    except OSError:
        print(f"Error opening file at path: {file_path}")
        exit()


# HERE STARTS THE PROGRAM
file_path = str(input("Insert your file path (.txt): "))  # Replace with your file path
sequences = read_fasta(file_path)
print(sequences)
# Ask the user how many sequences they want to align
num_sequences = int(input(f"\nHow many sequences do you want to align (1-{len(sequences)})? "))
if num_sequences > len(sequences) or num_sequences < 1:
    print("Invalid number of sequences.")
    exit()

selected_sequences = sequences[:num_sequences]
start_time = time.time()
aligned_sequences = multiple_sequence_alignment(selected_sequences)
end_time = time.time()
elapsed_time = end_time - start_time

# Measure memory usage
#mem_usage = memory_usage((multiple_sequence_alignment, (sequences,)), interval=0.1, timeout=1)

# Display the aligned sequences
print("\nAligned sequences:")
for seq in aligned_sequences:
    print(seq)

# Consensus sequences
a = str(input(f"\nWould you like a consensus sequence? (Y/N): ")).upper()

if a == "Y":
    dp = str(input("Are your sequences DNA or protein? (D/P): ")).upper()
    if dp == 'D':
        consensus_seq = get_consensus_dna(aligned_sequences)
    elif dp == 'P':
        consensus_seq = get_consensus_protein(aligned_sequences)
    else:
        print("Invalid input. Exiting.")
        exit()

    print("\nConsensus sequence:")
    print(consensus_seq)

print("Runtime: ",elapsed_time)
#print(f"Max Memory Usage: {max(mem_usage) if mem_usage else 0} MB")