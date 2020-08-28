#!/usr/bin/env nextflow

/*
========================================================================================
                         nextflow-lcep
========================================================================================
 https://github.com/mlf-core/nextflow-lcep
----------------------------------------------------------------------------------------
*/

normal_filtered_tpm_ch = Channel.fromPath(params.normal_tpm)
normal_replicates_ch = Channel.fromPath(params.normal_replicates)
cancer_filtered_tpm_ch = Channel.fromPath(params.cancer_tpm)
cancer_replicates_ch = Channel.fromPath(params.cancer_replicates)
pathways_ch = Channel.fromPath(params.pathways)

// TODO Refactor this and use tuples of files as input
process replicate_removal_normal {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file normal_tpm from normal_filtered_tpm_ch
    file normal_replicates from normal_replicates_ch
    file cancer_tpm from cancer_filtered_tpm_ch
    file cancer_replicates from cancer_replicates_ch

    output:
    file 'lung_normal_tpm_total_filtered_wo_rep.tsv' into normal_filtered_wo_replicates_ch
    file 'lung_cancer_tpm_total_filtered_wo_rep.tsv' into cancer_filtered_wo_replicates_ch

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
    file normal_lung_tpm from normal_filtered_wo_replicates_ch
    file cancer_lung_tpm from cancer_filtered_wo_replicates_ch
    file pathways from pathways_ch

    output:
    file 'human_pathways_tpm_lung_normal.tsv' into normal_pathways_tpm_ch
    file 'human_pathways_tpm_lung_cancer.tsv' into cancer_pathways_tpm_ch

    script:
    """
    reduce_to_kegg_pathways.py -n $normal_lung_tpm -c $cancer_lung_tpm -t lung -hp $pathways
    """
}

process generate_training_test_datasets {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file normal_pathways_tpm from normal_pathways_tpm_ch
    file cancer_pathways_tpm from cancer_pathways_tpm_ch

    output:
    file 'train.tsv' into training_data_ch
    file 'test.tsv' into test_data_ch

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
    file to_predict from test_data_ch

    output:
    file 'predictions.csv' into predicted_ch
    
    script:
    """
    lcep-package --input $to_predict --output predictions.csv
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