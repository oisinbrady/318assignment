from Bio import SeqIO
import pylab
import argparse


def plot_all(args: argparse.Namespace) -> None:
    """
    produce histograms of sequences lengths for each genome
    :param args: runtime arguments: -s/--s_min will exclude any sequence less than the user allocated number
    :return: plots graphs, written in visulization directory
    """

    in_dir = "/home/oisin/Documents/genomes/"
    fasta_file_names = {"RUG290.fa", "RUG404.fa",
                        "RUG700.fa", "RUG824.fa"}
    out_dir = "visualization/"
    for ff in fasta_file_names:
        genome_name = f"{ff[0:(ff.find('.'))]}"  # remove the suffix
        if args.s_min is not None:  # set the sequence length that we start recording records
            sl = args.s_min
        else:
            sl = 0
        record_sizes = [len(rec) for rec in SeqIO.parse(f"{in_dir}{ff}", "fasta")
                        if len(rec) > sl]
        print(f"\n'{genome_name}' contigs/scaffolds: {len(record_sizes)}")
        pylab.hist(record_sizes, bins=10)
        pylab.title(f"{len(record_sizes)} '{genome_name}' sequences\nLengths \
            {min(record_sizes)} to {max(record_sizes)}")
        pylab.xlabel("Sequence length (bp)")
        pylab.ylabel("Count")
        pylab.savefig(f"{out_dir}sequence_length_{genome_name}.png")
        print(f"\nresults written to: 'sequence_length_{genome_name}.png'")


def sequence_alignment() -> None:
    pass


def min_sequence_len(x) -> int:
    x = int(x)
    if x <= 0:
        raise argparse.ArgumentTypeError(
            "Minimum sequence length must be a positive integer")
    return x


def parse_argv() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Plot graphs for fasta genome files.')
    parser.add_argument("-s", "--s_min", type=min_sequence_len,
                        help="limit the minimum size \
                        of sequences to be plotted")
    return parser.parse_args()


def main():
    args = parse_argv()
    plot_all(args)


if __name__ == '__main__':
    main()
