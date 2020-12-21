import csv
import re
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm


IN_DIR = "/home/oisin/programs/cs318/318assignment/analysis/genome_annotation/function_prediction"
EGGNOG_FILES = ["/RUG290/job_MM_dv9sgs3x_annotations.tsv", "/RUG404/job_MM_roq44gaq_annotations.tsv"]


def parse_file():
    annotated_genes = dict()
    for file in EGGNOG_FILES:
        with open(f"{IN_DIR}{file}") as f:
            tsv_reader = csv.reader(f, delimiter='\t')
            annotations = [row for row in tsv_reader]

            print(annotations[3])  # field names
            print(annotations[4][20])  # COG functional category of first entry
            print("{:.8f}".format(float(annotations[4][2])))  # convert e-value to decimal

            annotated_genes[file] = annotations
    return annotated_genes


def normalised_freq_vector(genome_annotated_genes) -> dict:
    COG_GROUPS = ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Y', 'Z', 'A', 'B', 'J', 'K', 'L', 'C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q', 'R', 'S']
    # get the number of gene function predictions
    print(genome_annotated_genes[len(genome_annotated_genes)-3][0])
    q = re.findall(r'\d+', genome_annotated_genes[len(genome_annotated_genes)-3][0])
    total_queries = list(map(int, q))[0]
    print(total_queries)

    cog_count = dict()
    # count the occurances of each COG
    for i in range(4, total_queries):
        if genome_annotated_genes[i][20] in COG_GROUPS:
            cog_group = ''.join(genome_annotated_genes[i][20].split())
            if cog_group in cog_count:
                cog_count[cog_group] = cog_count[cog_group] + 1
            else:
                cog_count[cog_group] = 1

    # normalise the frequency value
    for cog_id in cog_count:
        cog_count[cog_id] = cog_count[cog_id] / total_queries

    return cog_count


def barchart(annotated_genes):
    COG_GROUPS = ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Y', 'Z', 'A', 'B', 'J', 'K', 'L', 'C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q', 'R', 'S']

    for genome in annotated_genes:
        cog_count = normalised_freq_vector(annotated_genes[genome])

        a = np.array([value for value in cog_count.values()])
        base_labels = [cog_id for cog_id in cog_count.keys()]  # x-axis labels
        y_pos = np.arange(len(base_labels))

        # color each bar uniquely
        cmap = cm.get_cmap('Accent')
        # an estimate of the range of frequency percentage for COGs
        norm = Normalize(vmin=0.0, vmax=0.25)

        # initialise the bar chart
        plt.bar(y_pos, a, align='center', alpha=0.5, color=cmap(norm(a)))
        plt.xticks(y_pos, base_labels)
        plt.ylabel("Normalised COG frequency")
        plt.ylim(0, 1.0)
        plt.xlabel("Function Class (COG ID)")
        plt.title(f"{genome[1:7]} COG Function Classification")

        out_dir = "/home/oisin/programs/cs318/318assignment/analysis/genome_annotation/function_prediction"
        plt.savefig(f"{out_dir}/{genome[1:7]}/{genome[1:7]}_barchart.png")
        plt.close()


def main():
    annotated_genes = parse_file()
    barchart(annotated_genes)


if __name__ == '__main__':
    main()