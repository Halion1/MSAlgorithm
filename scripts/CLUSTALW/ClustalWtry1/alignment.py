from lib.funciones import *
import json

# Archivo FASTA con las secuencias
fasta_file = "BVBRC_genome_sequence.fasta"

# Cantidad de secuencias para probar
sequence_counts = [5,7,9,11]
sequences_load = load_sequences(fasta_file)


alignment_times_acumulado = []
consensus_scores_acumulado  = []
longitudes_acumulado = {}
time_proceso_acumulado = []

datos = []

for x in range(3):
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
            datos.append({"Numero Secuencias":count, "Score": calculate_consensus_score(aligned_sequences), "Tiempo": alignment_time, "Longitudes": longitudes_alignment(aligned_sequences)})
        score = calculate_consensus_score(aligned_sequences)
        longitudes_acumulado[count]=longitudes_alignment(aligned_sequences)
        sequence_len.append(aligned_sequences)
        consensus_scores.append(score)

        

    alignment_times_acumulado.append(alignment_times)
    consensus_scores_acumulado.append(consensus_scores)
    end_time = time.time()

alignment_times_mean = np.mean(alignment_times_acumulado,axis=0)
alignment_times_std = np.std(alignment_times_acumulado,axis=0)

json_object = json.dumps({"Raw":datos}, indent=4)



with open("data.json", "w") as outfile:
    outfile.write(json_object)


# Output de los score y longitud
print("Estos son los scores",consensus_scores)
print("Estas son las longitudes de los genes",longitudes_acumulado)

# Crear el gráfico de barras
create_bar_chart(sequence_counts, alignment_times_mean, alignment_times_std, "Tiempo de Alineamiento", "Cantidad de Secuencias", "Tiempo (segundos)")
#create_bar_chart(sequence_counts, consensus_scores, "Score del Consenso", "Cantidad de Secuencias", "Score")

