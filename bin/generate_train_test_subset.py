#!/usr/bin/env python
import click
import datetime
import logging
import numpy as np
import sys
import time

from sklearn.model_selection import train_test_split

console = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


@click.command()
@click.option('-n', '--normal_samples', prompt='TPM table of normal samples',
              help='Location of the TPM table tsv file containing the normal tissue samples', required=True)
@click.option('-c', '--cancer_samples', prompt='TPM table of normal samples',
              help='Location of the TPM table tsv file containing the normal tissue samples', required=True)
def main(normal_samples, cancer_samples):
    # start timer
    start = time.time()

    # parse tpm values
    LOG.info("Parse TPM tables")
    transposed_counts_cancer, gene_names, sample_names_cancer = parse_counts_table(cancer_samples)
    transposed_counts_normal, gene_names, sample_names_normal = parse_counts_table(normal_samples)

    # transpose gene_names
    gene_names = [list(i) for i in zip(*gene_names)]

    # print number of samples to LOG
    LOG.info("Number cancer samples: " + str(len(transposed_counts_cancer)))
    LOG.info("Number normal samples: " + str(len(transposed_counts_normal)))

    # create list of true classes; 1 for cancer, 0 for normal
    LOG.info("Create list of true classes")
    true_classes = len(sample_names_cancer) * [1] + len(sample_names_normal) * [0]

    # merge normal and cancer info
    LOG.info("Merge normal and cancer samples to one dataset")
    all_tpm = transposed_counts_cancer + transposed_counts_normal
    all_sample_names = sample_names_cancer + sample_names_normal

    X_train, X_test, y_train, y_test, run_ids_train, run_ids_test = train_test_split(all_tpm, true_classes,
                                                                                     all_sample_names, test_size=0.25,
                                                                                     random_state=0,
                                                                                     stratify=true_classes)

    normal_filename = normal_samples.split("/")[-1].split(".")[0]
    cancer_filename = cancer_samples.split("/")[-1].split(".")[0]

    write_table_to_file(X_train, y_train, run_ids_train, gene_names, "train.tsv")
    write_table_to_file(X_test, y_test, run_ids_test, gene_names, "test.tsv")


def parse_counts_table(path):
    table = []
    gene_names = []
    with open(path, "r") as file:
        sample_names = next(file).split("\n")[0].split("\t")[2:]
        for line in file:
            splitted = line.split("\n")[0].split("\t")
            table.append([float(x) for x in splitted[2:]])
            gene_names.append(splitted[:2])

    table = [list(i) for i in zip(*table)]
    gene_names = [list(i) for i in zip(*gene_names)]

    return table, gene_names, sample_names


def write_table_to_file(x, y, sample_names, gene_names, outpath):
    x = [list(i) for i in zip(*x)]
    sample_names = [str(label) + "_" + sra_id for sra_id, label in zip(sample_names, y)]
    with open(outpath, "w") as file:
        file.write("Gene ID\tGene Name\t" + "\t".join(sample_names) + "\n")

        for gene_name, gene_tpms in zip(gene_names, x):
            gene_tpms = [str(x) for x in gene_tpms]
            file.write("\t".join(gene_name) + "\t")
            file.write("\t".join(gene_tpms) + "\n")


if __name__ == "__main__":
    sys.exit(main())
