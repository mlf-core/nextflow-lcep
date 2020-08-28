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
@click.option('-t', '--tpm', prompt='TPM table ',
              help='Location of the TPM table', required=True)
@click.option('-r', '--rep', prompt='replicate file',
              help='replicate file', required=True)
@click.option('-o', '--out', prompt='output file',
              help='output file', required=True)
def main(tpm, rep, out):
    # start timer
    start = time.time()

    # parse tpm values
    LOG.info("Parse TPM tables")
    transposed_counts, gene_names, sample_names = parse_counts_table(tpm)

    # print number of samples
    LOG.info("Number samples: " + str(len(transposed_counts)))


    LOG.info("Parse duplicate file")
    list_of_all_replicates, replicate_groups = parse_replicates(rep)

    replicates_median_tpms = []
    leading_ids = []
    count = 0
    for group in replicate_groups:
        tpm_values = []
        leading_id = ""
        for sra_id in group:
            try:
                idx = sample_names.index(sra_id)
                tpm_values.append(transposed_counts[idx])
                count += 1
                if leading_id == "":
                    leading_id = sra_id
            except:
                pass
        if not tpm_values == []:
            replicates_median_tpms.append(list(np.median(np.array(tpm_values), axis=0)))
            leading_ids.append(leading_id)

    output_counts_table = []
    ouput_sample_names = []
    for sample_tpms, sra_id in zip(transposed_counts, sample_names):
        if sra_id not in list_of_all_replicates:
            output_counts_table.append(sample_tpms)
            ouput_sample_names.append(sra_id)
    for sample_tpms, sra_id in zip(replicates_median_tpms, leading_ids):
        output_counts_table.append(sample_tpms)
        ouput_sample_names.append(sra_id)

    LOG.info("Number samples after replicate removal: " + str(len(output_counts_table)))

    output_counts_table = [list(i) for i in zip(*output_counts_table)]
    gene_names = [list(i) for i in zip(*gene_names)]
    
    with open(out, "w") as file:
        file.write("Gene ID\tGene Name\t" + "\t".join(ouput_sample_names) + "\n")
        for gene_name, gene_tpms in zip(gene_names, output_counts_table):
            gene_tpms = [str(x) for x in gene_tpms]
            file.write("\t".join(gene_name) + "\t")
            file.write("\t".join(gene_tpms) + "\n")


    # stop timer
    end = time.time()
    LOG.info("Script finished in " + str(round(end - start, 2)) + " sec")


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


def parse_replicates(path):
    list_of_all_replicates = list()
    replicate_groups = list()
    with open(path, "r") as file:
        for line in file:
            splitted = line.split("\n")[0].split("\t")
            replicate_groups.append(splitted)
            for sra_id in splitted:
                list_of_all_replicates.append(sra_id)

    return list_of_all_replicates, replicate_groups


if __name__ == "__main__":
    sys.exit(main())
