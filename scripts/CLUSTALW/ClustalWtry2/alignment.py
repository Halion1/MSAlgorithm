from lib.funciones import *
import json
import random
import multiprocessing

def processCadena(sequences_load,count):
    
    sequences_subset = sequences_load[:count]

    SeqIO.write(sequences_subset, "temp.fasta", "fasta")

    alignment_time = perform_alignment(sequences_subset)

    # Cargar el resultado del alineamiento
    with open("temp.aln", "r") as alignment_file:
        aligned_sequences = MultipleSeqAlignment(SeqIO.parse(alignment_file, "clustal"))
            
    score = calculate_consensus_score(aligned_sequences)
    longitudes = longitudes_alignment(aligned_sequences)
    # Guardado de resultados
    numero_aleatorio = random.randint(1000, 9999)
    nombre_archivo = f"Informes/Alineamientos/{count}--{numero_aleatorio}-{np.mean(longitudes)}.json"


    datos = {
        "Cadenas": longitudes, 
        "Alineamiento":{"Data": "","Tiempo": alignment_time,"Score": score, "Long": len(str(aligned_sequences[0].seq))} 
    }        
    json_object = json.dumps(datos, indent=4)

    with open(nombre_archivo, "w") as outfile:
        outfile.write(json_object)

if __name__ == '__main__':
    # Archivo FASTA con las secuencias
    fasta_file = "fastas/BVBRC_genome_sequence.fasta"

    # Cantidad de secuencias para probar
    sequence_counts = [5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    sequences_load = load_sequences(fasta_file)


    for numero in sequence_counts:
        procesos = []
        for x in range(10):
            proceso = multiprocessing.Process(target=processCadena, args=(sequences_load,numero,))
            procesos.append(proceso)
            proceso.start()

        # Esperar a que todos los procesos terminen
        for proceso in procesos:
            proceso.join()
        
        del procesos
'''
# Output de los score y longitud
print("Estos son los scores",consensus_scores)
print("Estas son las longitudes de los genes",longitudes_acumulado)

# Crear el gráfico de barras
create_bar_chart(sequence_counts, alignment_times_mean, alignment_times_std, "Tiempo de Alineamiento", "Cantidad de Secuencias", "Tiempo (segundos)")
#create_bar_chart(sequence_counts, consensus_scores, "Score del Consenso", "Cantidad de Secuencias", "Score")
'''
