from lib.funciones import *
import json
# Archivo FASTA con las secuencias
fasta_file = "secuencias (2).fasta"

#Carga Archivo
sequences_load = load_sequences(fasta_file)

#% de carga

fraccion_de_carga = [1,2,3]
cantidad_grupo = [2,3]

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
            temp_sequences_proces = []
            i0 = 0
            for i in range( len(Pseq.seq)//(frac+1),len(Pseq.seq),len(Pseq.seq)//(frac+1) ):
                temp_sequences_proces.append(Pseq[i0:i])
                i0 = i
            sequences_proces.append(temp_sequences_proces)
        
        aligned_sequences = []
        for idx in range(frac):
            # Separa N-ares para linealizar
            sequences_subset=[x[idx] for x in sequences_proces]
            SeqIO.write(sequences_subset, "temp.fasta", "fasta")

            Datos[can][frac]["Tiempo"].append(perform_alignment(sequences_subset))
            
            # Cargar el resultado del alineamiento
            
            with open("temp.aln", "r") as alignment_file:
                alignment = MultipleSeqAlignment(SeqIO.parse(alignment_file, "clustal"))
                Datos[can][frac]["Scores"].append(calculate_consensus_score(alignment))
                aligned_sequences.append(alignment)
        Datos[can][frac]["Seq"] = aligned_sequences
        Datos[can][frac]["Score"] = sum(Datos[can][frac]["Scores"])

print(Datos)




with open('resultados.txt', 'w') as f:
    for can in cantidad_grupo:
        for iterCam in range(can):
            #print(iterCam)
            for frac in fraccion_de_carga:
                for iterFrac in range(frac):
                    f. write(str(frac))
                    f. write(" ")
                    f.write(str(Datos[can][frac]["Seq"][iterFrac][iterCam].seq))
                f.write("\n")
            f.write("\n")
f.close




