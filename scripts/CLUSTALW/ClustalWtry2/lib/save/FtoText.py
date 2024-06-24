from Bio import SeqIO
# Nombre del archivo de texto y del archivo FASTA de salida
txt_file = "sequence_spike_protein.txt"
fasta_file = "sequences.fasta"

# Lee las secuencias del archivo de texto y crea un archivo FASTA
with open(txt_file, "r") as txt:
    records = []
    sequence_lines = txt.read().splitlines()
    for i, line in enumerate(sequence_lines):
        seq = line.strip()  # Elimina espacios en blanco al inicio y final
        record = f">Secuencia_{i+1}\n{seq}\n"  # Encabezado y secuencia
        records.append(record)

# Escribe las secuencias en el archivo FASTA
with open(fasta_file, "w") as fasta:
    fasta.writelines(records)

print("Archivo FASTA creado exitosamente.")

