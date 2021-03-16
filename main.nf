#!/usr/bin/env nextflow

/*
========================================================================================
                         nextflow-lcep
========================================================================================
 https://github.com/mlf-core/nextflow-lcep
----------------------------------------------------------------------------------------
*/

ch_liver_normal_tpm = Channel.fromPath(params.liver_normal_tpm)
ch_liver_cancer_tpm = Channel.fromPath(params.liver_cancer_tpm)
ch_pathways = Channel.fromPath(params.pathways)

process reduce_to_kegg_pathways {
    label 'with_cpus'
    publishDir "${params.outdir}/intermediate_results", mode: 'copy'

    input:
    file liver_normal_tpm from ch_liver_normal_tpm
    file liver_cancer_tpm from ch_liver_cancer_tpm
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
/*
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
*/