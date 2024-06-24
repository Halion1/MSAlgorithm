from lib.funciones import *
import json
# Archivo FASTA con las secuencias
fasta_file = "BVBRC_genome_sequence.fasta"

#Carga Archivo
sequences_load = load_sequences(fasta_file)

#% de carga

fraccion_de_carga = [100, 80, 70,60]
cantidad_grupo = [5,10,15,20]

Datos = {}

# Cargar Grupos
for can in cantidad_grupo:
    Datos[can] = {}
    for frac in fraccion_de_carga:
        Datos[can][frac] = {"Tiempo":[],"Seq":[],"Scores": [],"Score":0}
        # Franccionamiento de Senales
        sequences_pre_proces = sequences_load[:can]
        sequences_proces = []
        for Pseq in sequences_pre_proces:
            temp_sequences_proces = Pseq[:round(len(Pseq)*(frac/100))]
            
            sequences_proces.append(temp_sequences_proces)
        
        print(sequences_proces)
        SeqIO.write(sequences_proces, "temp.fasta", "fasta")

        Datos[can][frac]["Tiempo"].append(perform_alignment(sequences_proces))
        
        # Cargar el resultado del alineamiento
        
        with open("temp.aln", "r") as alignment_file:
            alignment = MultipleSeqAlignment(SeqIO.parse(alignment_file, "clustal"))
            Datos[can][frac]["Scores"].append(calculate_consensus_score(alignment))
            

        Datos[can][frac]["Seq"] = alignment
        Datos[can][frac]["Score"] = sum(Datos[can][frac]["Scores"])
print()
print(Datos)

with open('resultados.txt', 'w') as f:
    for can in cantidad_grupo:
        for iterCam in range(can):
            for frac in fraccion_de_carga:
                f. write(str(frac))
                f. write("\t")
                f.write(str(Datos[can][frac]["Seq"][iterCam].seq))
                f.write("\n")
            f.write("\n")
f.close
