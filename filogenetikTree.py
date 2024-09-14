from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# Membaca file alignment dalam format FASTA
alignment = AlignIO.read("AllignmentSpecies.aln-fasta", "fasta")

# Menggunakan DistanceCalculator untuk menghitung jarak (menggunakan model identitas)
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Membangun pohon dengan metode Neighbor-Joining (NJ)
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

# Menampilkan pohon filogenetik secara grafis
Phylo.draw(tree)

# Jika ingin menyimpan pohon dalam file Newick
Phylo.write(tree, "pohon_filo.newick", "newick")
