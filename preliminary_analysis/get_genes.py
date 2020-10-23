from Bio import SeqIO
import pyrodigal

# get all predicted and translated coding regions for a sequence
in_dir = "/home/oisin/programs/cs318/assignment/genomes/"
records = [rec for rec in SeqIO.parse(f"{in_dir}RUG290.fa", "fasta")]
p = pyrodigal.Pyrodigal(meta=True)
coding_regions = dict()
for record in records:
    proteins = list()
    for gene in p.find_genes(str(record.seq)):
        # TODO offer an option to gather the untranslated genes aswell
        proteins.append(gene.translate())
    coding_regions[record.id] = proteins
print(coding_regions)

