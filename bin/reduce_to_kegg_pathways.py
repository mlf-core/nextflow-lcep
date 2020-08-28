#!/usr/bin/env python
import click
import datetime
import logging
import numpy as np
import sys
import time

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
@click.option('-t', '--tissue', prompt='tissue name', help='tissue name', required=True)
@click.option('-hp', '--human_pathways', prompt='human pathways',
              help='all human pathways (KEGG)', required=True)
def main(normal_samples, cancer_samples, tissue, human_pathways):
    # start timer
    start = time.time()

    # parse tpm values
    LOG.info("Parse TPM tables")
    transposed_counts_cancer, gene_names, sample_names_cancer = parse_counts_table(cancer_samples)
    transposed_counts_normal, gene_names, sample_names_normal = parse_counts_table(normal_samples)

    # print number of samples to LOG
    LOG.info("Number cancer samples: " + str(len(transposed_counts_cancer)))
    LOG.info("Number normal samples: " + str(len(transposed_counts_normal)))

    LOG.info("Convert to np array")
    transposed_counts_cancer = np.array(transposed_counts_cancer)
    transposed_counts_normal = np.array(transposed_counts_normal)

    # parse all human pathway_genes
    human_pathway_genes = parse_genes_from_kegg(human_pathways)

    # Reduce genes to human pathways (KEGG)
    LOG.info("Feature reduction based on KEGG pathways")
    human_pathways_cancer_tpm, gene_names_human_pathways = reduce_genes(transposed_counts_cancer, gene_names,
                                                                        human_pathway_genes)
    human_pathways_normal_tpm, gene_names_human_pathways = reduce_genes(transposed_counts_normal, gene_names,
                                                                        human_pathway_genes)

    LOG.info("Number of genes after feature reduction all human pathways: " + str(len(human_pathways_cancer_tpm[0])))

    LOG.info("Write reduced TPM tables to file")
    outpath_human_pathway_cancer = "human_pathways_tpm_" + tissue + "_cancer.tsv"
    outpath_human_pathway_normal = "human_pathways_tpm_" + tissue + "_normal.tsv"

    generate_tpm_table(outpath_human_pathway_cancer, human_pathways_cancer_tpm, gene_names_human_pathways,
                       sample_names_cancer)
    generate_tpm_table(outpath_human_pathway_normal, human_pathways_normal_tpm, gene_names_human_pathways,
                       sample_names_normal)

    # stop timer
    end = time.time()
    LOG.info("Script finished in " + str(round(end - start, 2)) + " sec")


def reduce_genes(tpm_table, gene_names, pathway_genes):
    bool_mask = []
    for name in gene_names[0]:
        if name in pathway_genes:
            bool_mask.append(True)
        else:
            bool_mask.append(False)
    bool_mask = np.array(bool_mask)

    tpm_table = tpm_table[:, bool_mask]
    gene_name_copy = np.array(gene_names)

    gene_name_copy = gene_name_copy[:, bool_mask]

    return tpm_table, gene_name_copy


def parse_genes_from_kegg(pathway_file):
    pathway_name = []
    pathway_genes = []
    unique_genes = set()
    with open(pathway_file, "r") as file:
        for line in file:
            splitted = line.split("\n")[0].split("\t")
            if len(splitted) >= 3:
                pathway_name.append(splitted[0])
                pathway_genes.append(splitted[1:])
                for gene_id in splitted[1:]:
                    unique_genes.add(gene_id)
    unique_genes = list(unique_genes)
    return unique_genes


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


def generate_tpm_table(output_path, transposed_tpm, gene_names, sample_names):
    tpm_values = [list(i) for i in zip(*list(transposed_tpm))]
    gene_names = [list(i) for i in zip(*list(gene_names))]

    gene_names = list(gene_names)
    with open(output_path, "w") as file:
        file.write("Gene ID\tGene Name\t" + "\t".join(sample_names) + "\n")
        for gene_tpms, gene_name in zip(tpm_values, gene_names):

            gene_tpms = [str(x) for x in gene_tpms]
            file.write("\t".join(gene_name) + "\t" + "\t".join(gene_tpms) + "\n")


if __name__ == "__main__":
    sys.exit(main())
