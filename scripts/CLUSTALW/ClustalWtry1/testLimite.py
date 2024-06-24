from lib.funciones import *

# Archivo FASTA con las secuencias
fasta_file = "BVBRC_genome_sequence.fasta"

# Cantidad de secuencias para probar
sequence_counts = [2,3,4,5]
sequences_load = load_sequences(fasta_file)




alignment_times_acumulado = []
consensus_scores_acumulado  = []
longitudes_acumulado = {}
time_proceso_acumulado = []

for x in range(1):
    start_time = time.time()
    # Tiempos de alineamiento y scores
    alignment_times = []
    consensus_scores = []
    sequence_len  = []
    for count in sequence_counts:

        sequences_subset = sequences_load[:count]
            
        SeqIO.write(sequences_subset, "temp.fasta", "fasta")



        alignment_time = perform_alignment(sequences_subset)
        alignment_times.append(alignment_time)




        # Cargar el resultado del alineamiento
        aligned_sequences = []
        with open("temp.aln", "r") as alignment_file:
            alignment = MultipleSeqAlignment(SeqIO.parse(alignment_file, "clustal"))
            aligned_sequences = alignment

        score = calculate_consensus_score(aligned_sequences)
        longitudes_acumulado[count]=longitudes_alignment(aligned_sequences)
        sequence_len.append(aligned_sequences)
        consensus_scores.append(score)

        with open('Limite.txt', 'a') as f:
            f. write("Grupo de "+str(count)+"-> Puntaje: "+str(score) +", Longitud: "+str(longitudes_acumulado[count]) +", Tiempo: "+str(alignment_time)+"\n")
            f.close

    alignment_times_acumulado.append(alignment_times)
    consensus_scores_acumulado.append(consensus_scores)
    end_time = time.time()

    

