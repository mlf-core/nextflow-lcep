#!/usr/bin/env nextflow

/*
========================================================================================
                         nextflow-lcep
========================================================================================
 https://github.com/mlf-core/nextflow-lcep
----------------------------------------------------------------------------------------
*/

ch_normal_filtered_tpm = Channel.fromPath(params.normal_tpm)
ch_normal_replicates = Channel.fromPath(params.normal_replicates)
ch_cancer_filtered_tpm = Channel.fromPath(params.cancer_tpm)
ch_cancer_replicates = Channel.fromPath(params.cancer_replicates)
ch_pathways = Channel.fromPath(params.pathways)

// TODO Refactor this and use tuples of files as input
process replicate_removal_normal {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file normal_tpm from ch_normal_filtered_tpm
    file normal_replicates from ch_normal_replicates
    file cancer_tpm from ch_cancer_filtered_tpm
    file cancer_replicates from ch_cancer_replicates

    output:
    file 'lung_normal_tpm_total_filtered_wo_rep.tsv' into ch_normal_filtered_wo_replicates
    file 'lung_cancer_tpm_total_filtered_wo_rep.tsv' into ch_cancer_filtered_wo_replicates

    script:
    """
    replicate_removal.py -t $normal_tpm -r $normal_replicates -o lung_normal_tpm_total_filtered_wo_rep.tsv
    replicate_removal.py -t $cancer_tpm -r $cancer_replicates -o lung_cancer_tpm_total_filtered_wo_rep.tsv
    """
}

process reduce_to_kegg_pathways {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file normal_lung_tpm from ch_normal_filtered_wo_replicates
    file cancer_lung_tpm from ch_cancer_filtered_wo_replicates
    file pathways from ch_pathways

    output:
    file 'human_pathways_tpm_lung_normal.tsv' into ch_normal_pathways_tpm
    file 'human_pathways_tpm_lung_cancer.tsv' into ch_cancer_pathways_tpm

    script:
    """
    reduce_to_kegg_pathways.py -n $normal_lung_tpm -c $cancer_lung_tpm -t lung -hp $pathways
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
    lcep-package --input $to_predict --output predictions.csv --cuda
    """
}

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
