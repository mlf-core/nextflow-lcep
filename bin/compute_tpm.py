#!/usr/bin/env python
import click
import sys
import time
import logging
import statistics
from multiprocessing import Pool
from tqdm import tqdm

console = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("TPM calculator")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inputfile', prompt='input file',
              help='Input file; Count table in the format of the recount project',
              required=True)
@click.option('-o', '--outputfile', prompt='output file',
              help='Output file; Table containing TPM values',
              required=True)
@click.option('-p', '--processes',
              help='number of processes to use for the tpm computation',
              required=False, default=4)
@click.option('-g', '--genelength_file', prompt='gene length file',
              help='gene length file produced by gtftools.py',
              required=True)
def main(inputfile, outputfile, processes, genelength_file):
    # start timer
    start = time.time()

    LOG.info("Parsing counts table")
    alt_table = parse_data(inputfile)

    # parse gene length file
    global gene_lengths
    gene_lengths = parse_gene_length_file(genelength_file)

    # number of samples
    num_samples = len(alt_table[0]) - 2  # two gene name/ID columns

    # transpose array
    transposed_alt_table = [list(i) for i in zip(*alt_table)]
    global gene_names
    gene_names = transposed_alt_table[:2]
    transposed_counts = transposed_alt_table[2:]

    LOG.info("TPM computation")
    #pool = Pool(processes=processes)
    #tpm_table = pool.map(compute_single_tpm, transposed_counts)
    tpm_table = []
    for gene_counts in tqdm(transposed_counts):
        tpm_table.append(compute_single_tpm(gene_counts))

    # transpose back
    tpm_table = [list(i) for i in zip(*tpm_table)]
    gene_names = [list(i) for i in zip(*gene_names)]

    # write to file
    LOG.info("Generating output file")
    write_to_file(tpm_table, gene_names, outputfile)

    # stop timer
    end = time.time()
    LOG.info("Process finished in " + str(round(end - start, 2)) + " sec for " + str(num_samples) + " samples")

    return 0


def parse_data(inputfile):
    """
    :param inputfile: path to file
    :return: 2D array counts and gene names; first two cols gene names, test counts
    """
    table = list()
    with open(inputfile, "r") as file:
        for line in file:
            splitted = line.split("\n")[0].split("\t")
            table.append(splitted)
    return table


def parse_gene_length_file(path):
    """
    :param path: path to GTF file
    :return: dict (gene_id: median transcript length)
    """
    gene_length_dict = {}
    with open(path, "r") as file:
        next(file)
        for line in file:
            splitted_line = line.split("\n")[0].split("\t")
            ensembl_id = splitted_line[0]
            median_length = splitted_line[2]
            gene_length_dict[ensembl_id] = int(median_length)
    return gene_length_dict


def compute_single_tpm(sample):
    """
    :param sample: counts of single samples
    :return: tpm values
    """
    sample_sum = 0
    for i in range(1, len(sample)):
        sample_sum += int(float(sample[i])) / gene_lengths[gene_names[0][i]]
    for j in range(1, len(sample)):
        if not sample[j] == "0":
            sample[j] = str(((int(float(sample[j])) / gene_lengths[gene_names[0][j]]) * 1000000) / sample_sum)
    return sample


def write_to_file(tpm_values, gene_names, outputfile):
    """
    :param tpm_values: tpm table
    :param gene_names: two columns with gene ID and gene name
    :param outputfile: path for output file
    :return: None
    """
    with open(outputfile, "w") as file:
        for i in range(len(tpm_values)):
            file.write("\t".join(gene_names[i]) + "\t" + "\t".join(tpm_values[i]) + "\n")


if __name__ == "__main__":
    sys.exit(main())
