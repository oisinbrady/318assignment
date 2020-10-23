from Bio import SeqIO


def fasta_to_dict():
    # convert all 4 genome sequences into Biopython sequence data types
    file_dir = "/home/oisin/programs/cs318/assignment/genomes/"
    fasta_file_names = {"RUG290.fa", "RUG404.fa",
                        "RUG700.fa", "RUG824.fa"}
    all_genomes = dict()
    for file in fasta_file_names:
        with open(f"{file_dir}{file}", "r") as ff:
            sequences = list(SeqIO.parse(ff, "fasta"))
            all_genomes[f"{file[0:(file.find('.'))]}"] = sequences
    return all_genomes


def count_n_bases(all_genomes: dict):
    for genome in all_genomes:
        a, g, c, t = 0, 0, 0, 0
        for record in all_genomes[genome]:
            a += record.seq.count('A')
            g += record.seq.count('G')
            c += record.seq.count('C')
            t += record.seq.count('T')
        print(f"\nGenome:'{genome}':\
            \nNucleotide count:\
            \nA:{a}, G:{g}, C:{c}, T:{t}\n")


def main():
    all_genomes = fasta_to_dict()
    count_n_bases(all_genomes)


if __name__ == '__main__':
    main()
