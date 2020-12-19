from Bio import SeqIO
from itertools import permutations
from math import sqrt
from matplotlib import pyplot as plt
import numpy as np
import sys


def fasta_to_dict(genome_file=None):
    if genome_file == None:
        # convert all 4 genome sequences into Biopython sequence data types
        file_dir = "/home/oisin/programs/cs318/318assignment/genomes/"
        fasta_file_names = {"RUG290.fa", "RUG404.fa",
                            "RUG700.fa", "RUG824.fa"}
        all_genomes = dict()
        for file in fasta_file_names:
            with open(f"{file_dir}{file}", "r") as ff:
                sequences = list(SeqIO.parse(ff, "fasta"))
                all_genomes[f"{file[0:(file.find('.'))]}"] = sequences
        return all_genomes
    else:
        genome = dict()
        with open(f"{file_dir}{file}", "r") as ff:
            sequences = list(SeqIO.parse(ff, "fasta"))
            genome[f"{file[0:(file.find('.'))]}"] = sequences
        return genome


def count_n_bases(all_genomes: dict, k: int):
    k_mers = list(permutations(['A', 'G', 'C', 'T', 'N'], k))
    all_kmer_vectors = dict()
    for genome in all_genomes:
        # init dictionary of k-mer counts
        k_mer_count = dict()
        for kmer in k_mers:
            kmer_id = "".join(kmer)
            k_mer_count[kmer_id] = 0

        # count each k-mer permutation
        for record in all_genomes[genome]:
            for kmer in k_mers:
                kmer_id = "".join(kmer)
                k_mer_count[f"{kmer_id}"] = k_mer_count.get(kmer_id) + record.seq.count(f"{kmer_id}")

        print(f"\nGenome:'{genome}':\nk-mer count:{k_mer_count}")

        kmer_vector = normalised_k_mer_freq_vector(k_mer_count, k)
        all_kmer_vectors[f"{genome}"] = kmer_vector
        print(f"Normalised {k}-mer frequency vector:\n{kmer_vector}\n")

    return all_kmer_vectors


def normalised_k_mer_freq_vector(k_mer_count: dict, k: int):
    norm_k_mer_freq_vector = dict()
    total_bases = sum(k_mer_count.values())
    for kmer in k_mer_count:
        norm_k_mer_freq_vector[kmer] = k_mer_count.get(kmer)/total_bases
    print(f"Total Nucleotides: {total_bases}")
    return norm_k_mer_freq_vector


def euclidean_distance(kmer_vectors: dict, k: int):
    """
    calculate euclidean distance between all genomes
    N.B: Ns (unidentified bases) are exluded in this calculation
    """
    print("Euclidean Distances:\n")
    complete_v = list()  # vectors already compared
    for cur_v in kmer_vectors:
        for next_v in kmer_vectors:
            if next_v != cur_v:  # compare with all other genome vectors
                if next_v not in complete_v:
                    # get normalised value of each kmer of the next genome
                    kmer_delta_vals = list()
                    for item in kmer_vectors[next_v]:
                        """
                        calculate the difference in value of the same kmer
                        between the current genome and the next genome
                        """
                        cur_kmer_val = kmer_vectors[cur_v].get(item)
                        next_kmer_val = kmer_vectors[next_v].get(item)
                        # delta value is squared for euclidean distance calc
                        kmer_delta = (cur_kmer_val - next_kmer_val) ** 2
                        kmer_delta_vals.append(kmer_delta)
                    eucl_d = sqrt(sum(kmer_delta_vals))
                    print(f"d({cur_v}, {next_v})) = {eucl_d}")

        # memoize vectors already compared to avoid redundanies in iteration
        complete_v.append(cur_v)


def histogram(kmer_vectors: dict) -> None:
    for genome_name in kmer_vectors:
        cur_v = kmer_vectors[genome_name]
        dataset = list()
        for item in cur_v:
            dataset.append(cur_v.get(item))
        a = np.array(dataset)
        base_labels = [item for item in cur_v]
        y_pos = np.arange(len(base_labels))

        plt.bar(y_pos, a, align='center', alpha=0.5)
        plt.xticks(y_pos, base_labels)
        plt.ylabel("normalised frequency")
        plt.xlabel("k-mer")
        plt.title(genome_name)

        out_dir = "/home/oisin/programs/cs318/318assignment/analysis/kmer_analysis/histograms"
        plt.savefig(f"{out_dir}/{genome_name}_hist.png")
        plt.close()


def main():
    if len(sys.argv) > 2:
        # i.e., if a file is also passed
        k = int(sys.argv[1])
        all_genomes = fasta_to_dict(sys.argv[2])
    elif len(sys.argv) > 1:
        k = int(sys.argv[1])
        all_genomes = fasta_to_dict()
    else:
        k = 1  # default k-mer size

    all_genome_kmer_vectors = count_n_bases(all_genomes, k)
    euclidean_distance(all_genome_kmer_vectors, k)
    histogram(all_genome_kmer_vectors)


if __name__ == '__main__':
    main()
