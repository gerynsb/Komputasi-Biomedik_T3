# Fungsi untuk membaca file FASTA untuk mengetahui DNA pada Anoa
def read_genome_fasta(FullGenome):
    with open(FullGenome, 'r') as file:
        fasta_sequence = ''.join([line.strip() for line in file if not line.startswith('>')])
    return fasta_sequence

# Tambahkan argumen ke fungsi untuk membaca file DNA mitokondria spesifik pada anoa yang hanya sebagian DNA saja
def read_dna_fasta(anoaMithocondriaDNA):
    with open(anoaMithocondriaDNA, 'r') as file:
        fasta_sequence = ''.join([line.strip() for line in file if not line.startswith('>'  )])
    return fasta_sequence

# Membaca urutan genom keseluruhan dan DNA mitokondria
genome_sequence = read_genome_fasta('FullGenome.fna')
print('Genome keseluruhan dari anoa : ', genome_sequence)

dna_sequence = read_dna_fasta('anoaMithocondriaDNA.fasta')
print('\nDNA mitokondria dari anoa : ', dna_sequence)

# Fungsi untuk menerjemahkan DNA menjadi protein (asam amino)
def translate_dna_to_protein(sequence):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    protein_sequence = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            protein_sequence += codon_table.get(codon, '?')
    return protein_sequence

# Menerjemahkan genome dan DNA mitokondria menjadi protein
proteinGenome = translate_dna_to_protein(genome_sequence)
proteinDNA = translate_dna_to_protein(dna_sequence)

print("\nHasil asam amino genome: ", proteinGenome)
print("\nHasil asam amino DNA mitokondria: ", proteinDNA)


