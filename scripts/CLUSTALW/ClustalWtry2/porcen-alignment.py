from lib.funciones import *
import json
import random
import multiprocessing

def processCadena(sequences_load,count,frac):
    
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
    nombre_archivo = f"Informes/Porcentual/{count}--{numero_aleatorio}-{np.mean(longitudes)}.json"


    datos = {
        "Cadenas": longitudes, 
        "Alineamiento":{"Data": "","Tiempo": alignment_time,"Score": score, "Long": len(str(aligned_sequences[0].seq)), "Frac": frac} 
    }        
    json_object = json.dumps(datos, indent=4)

    with open(nombre_archivo, "w") as outfile:
        outfile.write(json_object)

def procentajeCadena(sequences_pre_proces,frac):
    sequences_proces = []
    for Pseq in sequences_pre_proces:
        temp_sequences_proces = Pseq[:round(len(Pseq)*(frac/100))]
        
        sequences_proces.append(temp_sequences_proces)
    return sequences_proces

if __name__ == '__main__':
    # Archivo FASTA con las secuencias
    fasta_file = "fastas/BVBRC_genome_sequence.fasta"

    # Cantidad de secuencias para probar
    sequence_counts = [19]
    sequences_porcentaje = [90,80,70,60,50,40,30,20]
    sequences_load_row = load_sequences(fasta_file)

    for procentaje in sequences_porcentaje:
        sequences_load = procentajeCadena(sequences_load_row,procentaje)
        for numero in sequence_counts:
            procesos = []
            for x in range(3):
                proceso = multiprocessing.Process(target=processCadena, args=(sequences_load,numero,procentaje,))
                procesos.append(proceso)
                proceso.start()

            # Esperar a que todos los procesos terminen
            for proceso in procesos:
                proceso.join()
            
            del procesos
    