#!/usr/bin/env nextflow
/*
========================================================================================
nextflow-lcep
========================================================================================
https://github.com/mlf-core/nextflow-lcep
----------------------------------------------------------------------------------------
*/
params.gtex_counts_link = "http://duffel.rail.bio/recount/v2/SRP012682/counts_gene.tsv.gz"
params.gtex_metadata_link = "http://duffel.rail.bio/recount/SRP012682/SRP012682.tsv"
params.tcga_counts_link = "http://duffel.rail.bio/recount/v2/TCGA/counts_gene.tsv.gz"
params.tcga_metadata_link = "http://duffel.rail.bio/recount/TCGA/TCGA.tsv"
params.gencode_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"

ch_gtex_counts_link = Channel.value(params.gtex_counts_link)
ch_gtex_metadata_link = Channel.value(params.gtex_metadata_link)
ch_tcga_counts_link = Channel.value(params.tcga_counts_link)
ch_tcga_metadata_link = Channel.value(params.tcga_metadata_link)
ch_gtf_file = Channel.value(params.gencode_link)
ch_pathways = Channel.fromPath(params.pathways)

process filter_gtex {
    label 'with_cpus'

    input:
    val gtex_metadata_file from ch_gtex_metadata_link
    val gtex_counts_file from ch_gtex_counts_link

    output:
    file 'filtered_gtex_counts.tsv' into ch_filtered_gtex_counts

    script:
    """
    wget $gtex_counts_file
    wget $gtex_metadata_file
    filter_gtex_counts.py -c counts_gene.tsv.gz -m SRP012682.tsv
    """
}

process filter_tcga {
    label 'with_cpus'

    input:
    val tcga_metadata_file from ch_tcga_metadata_link
    val tcga_counts_file from ch_tcga_counts_link

    output:
    file 'filtered_tcga_counts.tsv' into ch_filtered_tcga_counts

    script:
    """
    wget $tcga_counts_file
    wget $tcga_metadata_file
    filter_tcga_counts.py -c counts_gene.tsv.gz -m TCGA.tsv
    """
}

process sort_counts {
    label 'with_cpus'

    input:
    file filtered_tcga_counts from ch_filtered_tcga_counts
    file filtered_gtex_counts from ch_filtered_gtex_counts

    output:
    file 'healthy_counts.tsv' into ch_healthy_counts
    file 'cancer_counts.tsv' into ch_cancer_counts

    script:
    """
    sort_counts_types.py -g $filtered_gtex_counts -t $filtered_tcga_counts
    """
}

process get_gene_lengths {
    label 'with_cpus'

    input:
    val gtf from ch_gtf_file

    output:
    file 'gencode.v25.annotation.gtf.genelength' into ch_genelengths

    beforeScript 'chmod o+rw .'
    script:
    """
    wget $gtf
    gunzip gencode.v25.annotation.gtf.gz
    gtftools.py -l gencode.v25.annotation.gtf.genelength -c 1-22,X,Y,MT gencode.v25.annotation.gtf
    """
}

process tpm_conversion {
    label 'with_cpus'
    publishDir "${params.outdir}/tpm_conversion", mode: 'copy'
    cpus = 8
    memory '16 GB'

    input:
    file healthy_counts from ch_healthy_counts
    file cancer_counts from ch_cancer_counts
    file genelengths from ch_genelengths

    output:
    file 'healthy_tpm.tsv' into ch_healthy_tpm
    file 'cancer_tpm.tsv' into ch_cancer_tpm

    script:
    """
    compute_tpm.py -i $cancer_counts -o cancer_tpm.tsv -p 8 -g $genelengths
    compute_tpm.py -i $healthy_counts -o healthy_tpm.tsv -p 8 -g $genelengths
    """
}

process reduce_to_kegg_pathways {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file liver_normal_tpm from ch_healthy_tpm
    file liver_cancer_tpm from ch_cancer_tpm
    file pathways from ch_pathways

    output:
    file 'tpm_liver_cancerpathways_normal.tsv' into ch_normal_pathways_tpm
    file 'tpm_liver_cancerpathways_cancer.tsv' into ch_cancer_pathways_tpm

    script:
    """
    reduce_to_kegg_pathways.py -n $liver_normal_tpm -c $liver_cancer_tpm -t liver_cancerpathways -hp $pathways -o .
    """
}

process generate_training_test_datasets {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file normal_pathways_tpm from ch_normal_pathways_tpm
    file cancer_pathways_tpm from ch_cancer_pathways_tpm

    output:
    file 'train.tsv' into ch_training_data
    file 'test.tsv' into ch_test_data

    script:
    """
    generate_train_test_subset.py -n $normal_pathways_tpm -c $cancer_pathways_tpm
    """
}

process predict_lcep {
    echo true
    label 'with_all_gpus'
    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    file to_predict from ch_test_data

    output:
    file 'predictions.csv' into ch_predicted

    script:
    """
    lcep-package --input $to_predict --output predictions.csv
    """
}
/*
process run_system_intelligence {
    publishDir "${params.outdir}/results", mode: 'copy'
    label 'with_all_gpus'

    output:
    file 'system_intelligence.html'
    file 'system_intelligence.json'

    script:
    """
    system-intelligence all --output_format json --generate_html_table --output system_intelligence.json
    """
}
*/
