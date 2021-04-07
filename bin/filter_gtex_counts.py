#!/usr/bin/env python
import click
import datetime
import logging
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

    # scan for liver sample IDs
    LOG.info("Parse metadata file")
    liver_ids = parse_metadata(metadata)

    LOG.info("filter liver counts")
    filtered_counts, gene_names, ids_liver_header = filter_counts(liver_ids, counts)

    LOG.info("Generate filtered output count table")
    write_count_table(filtered_counts, gene_names, ids_liver_header)

    # stop timer
    end = time.time()
    LOG.info("Script finished in " + str(round(end - start, 2)) + " sec")


def parse_metadata(metadata):
    """
    Searches the file for rows which have "Liver" specified as tissue

    :param metadata: path to metadata file (phenotype file of the recount2 homepage)
    :return: list of SRA ids of liver samples
    """
    liver_ids = []
    with open(metadata, "r") as file:
        header = next(file)
        for line in file:
            try:
                splitted_line = line.split("\n")[0].split("\t")
                run_id = splitted_line[3]
                tissue = splitted_line[25]
                if tissue == "Liver":
                    liver_ids.append(run_id)
            except:  # in case of inconsistent line formatting
                pass
    return liver_ids


def filter_counts(liver_ids, counts):
    """
    Filters columns which have the SRR ID of the given list

    :param liver_ids: list of SRA ids of liver samples
    :param counts: path to counts table
    :return:
    """
    filtered_counts = []
    gene_names = []
    with gzip.open(counts, 'r') as file:
        header = next(file).decode("utf-8")
        header_splitted = header.split("\n")[0].split("\t")
        ids_file = header_splitted[:-1]

        bool_array_file = [True if my_id in liver_ids else False for my_id in ids_file]
        ids_liver_header = list(compress(ids_file, bool_array_file))

        for line in file:
            splitted = line.decode("utf-8").split("\n")[0].split("\t")
            count_line = splitted[:-1]
            filtered_count_line = list(compress(count_line, bool_array_file))
            gene_name = splitted[-1]
            gene_names.append(gene_name)
            filtered_counts.append(filtered_count_line)

    ids_liver_header = ["0_" + x for x in ids_liver_header]
    return filtered_counts, gene_names, ids_liver_header


def write_count_table(filtered_counts, gene_names, ids_liver_header):
    """
    Writes information to file

    :param filtered_counts: for each gene, a list of gene count of each file (nested)
    :param gene_names: list of gene names
    :param ids_liver_header: SRR ids (appended with "_0" indicating healthy tissue
    :return:
    """
    with open("filtered_gtex_counts.tsv", "w") as file:
        file.write("gene_id\tgene_name\t" + "\t".join(ids_liver_header) + "\n")
        for gene_name, counts_line in zip(gene_names, filtered_counts):
            file.write(gene_name + "\t" + "\t" + "\t".join(counts_line) + "\n")


if __name__ == "__main__":
    sys.exit(main())
