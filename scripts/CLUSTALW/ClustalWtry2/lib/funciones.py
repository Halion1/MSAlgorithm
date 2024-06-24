import time
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

#Incio codigo

# Cargar las secuencias desde un archivo FASTA
def load_sequences(file_name):
    sequences = []
    with open(file_name, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            #print(record)
            sequences.append(record)
    return sequences

# Realizar el alineamiento utilizando ClustalW y medir el tiempo
def perform_alignment(sequences):
    start_time = time.time()
    #print("Antes de w")
    alignment_clustalw = ClustalwCommandline("clustalw", infile="temp.fasta")
    print(alignment_clustalw())
    #print("Despues de w")
    end_time = time.time()
    alignment_time = end_time - start_time
    return alignment_time

def longitudes_alignment(alignment):
    contador = []
    sequences = [len(str(record.seq))-str(record.seq).count('-') for record in alignment]
    return sequences


# Calcular el score del consenso del alineamiento
def calculate_consensus_score(alignment):
    contador = []
    sequences = [str(record.seq) for record in alignment]
    sequences_unit = []
    for row in sequences:
        sequences_unit.append([x for x in row])
    sequences_transpo=np.transpose(sequences_unit)
    
    for row in sequences_transpo:
        unique = np.unique(row)
        
        if( "-" in unique ):
            contador.append(-1)
        elif( len(unique)>1 ):
            contador.append(0)
        else:
            contador.append(1)
    
    return sum(contador)/max(longitudes_alignment(alignment))

# Crear un gráfico de barras para diferentes cantidades de secuencias
def create_bar_chart(x_labels, y_values , std, title, xlabel, ylabel):
    plt.bar(x_labels, y_values,yerr=std)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()