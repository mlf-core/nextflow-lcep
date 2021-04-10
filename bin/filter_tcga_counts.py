#!/usr/bin/env python
import click
import datetime
import logging
import numpy as np
import sys
import time
import gzip
from itertools import compress

console = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


@click.command()
@click.option('-m', '--metadata', prompt='GTEx metadata file',
              help='GTEx metadata (phenotype) file of the recount2 project', required=True)
@click.option('-c', '--counts', prompt='GTEx count table (tsv.gz)',
              help='GTEx count table of the recount2 project (tsv.gz)', required=True)
def main(metadata, counts):
    # start timer
    start = time.time()

    # scan for TCGA-LIHC sample IDs
    LOG.info("Parse metadata file")
    id_to_classes_recount = parse_metadata(metadata)
    tcga_lihc_ids = [x.upper() for x in id_to_classes_recount.keys()] #

    LOG.info("filter TCGA counts")
    filtered_counts, gene_names, ids_tcga_header = filter_counts(tcga_lihc_ids, counts)

    # add 0/1 class to sample ID
    ids_tcga_header_class = [str(id_to_classes_recount[x]) + "_" + x for x in ids_tcga_header]

    LOG.info("Generate filtered output count table")
    write_count_table(filtered_counts, gene_names, ids_tcga_header_class)

    # stop timer
    end = time.time()
    LOG.info("Script finished in " + str(round(end - start, 2)) + " sec")


def parse_metadata(metadata):
    """
    Searches the file for rows which have "TCGA-LIHC" specified as project and identifies tumor/healthy tissue samples

    :param metadata: path to metadata file (phenotype file of the recount2 homepage)
    :return: dictionary with IDs as key and 0/1 as class (healthy = 0, tumor = 1)
    """
    id_to_classes_recount = {}
    with open(metadata, "r") as file:
        header = next(file)
        for line in file:
            try:
                splitted_line = line.split("\n")[0].split("\t")
                file_id = splitted_line[22]
                project = splitted_line[77]
                sample_type = splitted_line[107]
                if project == "TCGA-LIHC":
                    if sample_type == 'Primary Tumor':
                        id_to_classes_recount[file_id] = 1
                    elif sample_type == 'Solid Tissue Normal':
                        id_to_classes_recount[file_id] = 0
                    elif sample_type == 'Recurrent Tumor':
                        id_to_classes_recount[file_id] = 1
                    else:
                        print(sample_type)
            except:
                pass
    return id_to_classes_recount


def filter_counts(tcga_lihc_ids, counts):
    filtered_counts = []
    gene_names = []
    with gzip.open(counts, 'r') as file:
        header = next(file).decode("utf-8")
        header_splitted = header.split("\n")[0].split("\t")
        ids_file = header_splitted[:-1]

        bool_array_file = [True if my_id in tcga_lihc_ids else False for my_id in ids_file]
        ids_tcga_header = list(compress(ids_file, bool_array_file))
        ids_tcga_header = [x.lower() for x in ids_tcga_header]

        for line in file:
            splitted = line.decode("utf-8").split("\n")[0].split("\t")
            count_line = splitted[:-1]
            filtered_count_line = list(compress(count_line, bool_array_file))
            gene_name = splitted[-1]
            gene_names.append(gene_name)
            filtered_counts.append(filtered_count_line)

    return filtered_counts, gene_names, ids_tcga_header


def write_count_table(filtered_counts, gene_names, ids_liver_header):
    """
    Writes information to file

    :param filtered_counts: for each gene, a list of gene count of each file (nested)
    :param gene_names: list of gene names
    :param ids_liver_header: TCGA sample IDs  (appended with "0_" or "1_" indicating healthy/tumor tissue )
    :return:
    """
    with open("filtered_tcga_counts.tsv", "w") as file:
        file.write("gene_id\tgene_name\t" + "\t".join(ids_liver_header) + "\n")
        for gene_name, counts_line in zip(gene_names, filtered_counts):
            file.write(gene_name + "\t" + "\t" + "\t".join(counts_line) + "\n")


if __name__ == "__main__":
    sys.exit(main())
