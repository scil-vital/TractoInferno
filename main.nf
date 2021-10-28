#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas":"$params.atlas",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

// Hardcoded Recobundle segmentation parameters
multi_parameters="18"
minimal_vote_ratio="0.5"
wb_clustering_thr="15 12"
seeds="0"
outlier_alpha="0.5"

// Atlas config
atlas_directory="$params.atlas/rbx_atlas"
atlas_anat="$params.atlas/mni_masked.nii.gz"
atlas_config="$params.atlas/config.json"

log.info "SCIL TractoInferno evaluation pipeline"
log.info "=========================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Atlas]"
log.info "Atlas Config: $atlas_config"
log.info "Atlas Anat: $atlas_anat"
log.info "Atlas Directory: $atlas_directory"
log.info ""
log.info "[Recobundles options]"
log.info "Multi-Parameters Executions: $multi_parameters"
log.info "Minimal Vote Percentage: $minimal_vote_ratio"
log.info "Whole Brain Clustering Threshold: $wb_clustering_thr"
log.info "Random Seeds: $seeds"
log.info "Outlier Removal Alpha: $outlier_alpha"
log.info ""
log.info ""

log.info "Input: $params.input"
log.info "Reference: $params.reference"
root = file(params.input)
ref = file(params.reference)


/* Watch out, files are ordered alphabetically in channel */
tractogram_for_recognition = Channel
    .fromFilePairs("$root/**/{*.*,}",
                    size:-1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromPath("$ref/**/*fa.nii.gz",
                    maxDepth:2)
    .map{[it.parent.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_bundles}

ref_bundles = Channel
    .fromFilePairs("$ref/**/{*.trk,}",
                    size:-1,
                    maxDepth:2) {it.parent.parent.name}

Channel.fromPath("$atlas_anat")
    .into{atlas_anat;atlas_anat_for_average}
atlas_config = Channel.fromPath("$atlas_config")
atlas_directory = Channel.fromPath("$atlas_directory")


workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_average
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__native_anat.nii.gz"
    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
    cp ${native_anat} ${sid}__native_anat.nii.gz
    """
}


tractogram_for_recognition
    .join(anat_for_reference_bundles)
    .join(transformation_for_recognition)
    .combine(atlas_config)
    .combine(atlas_directory)
    .set{tractogram_and_transformation}

process Recognize_Bundles {
    cpus params.rbx_processes
    memory params.rbx_memory_limit

    input:
    set sid, file(tractograms), file(reference), file(transfo), file(config), file(directory) from tractogram_and_transformation
    output:
    set sid, "*.trk" optional true into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"
    script:
    """
    if [ `echo $tractograms | wc -w` -gt 1 ]; then
        scil_streamlines_math.py lazy_concatenate $tractograms tracking_concat.trk
    else
        mv $tractograms tracking_concat.trk
    fi
    scil_remove_invalid_streamlines.py tracking_concat.trk tractogram_ic.trk --reference ${reference} --remove_single_point --remove_overlapping_points
    mkdir tmp/
    scil_recognize_multi_bundles.py tractogram_ic.trk ${config} ${directory}/*/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --multi_parameters $multi_parameters --minimal_vote_ratio $minimal_vote_ratio \
        --tractogram_clustering_thr $wb_clustering_thr --seeds $seeds --processes $params.rbx_processes
    rm tractogram_ic.trk tracking_concat.trk
    mv tmp/* ./
    """
}

bundles_for_cleaning
    .transpose()
    .set{all_bundles_for_cleaning}
process Clean_Bundles {
    input:
    set sid, file(bundle) from all_bundles_for_cleaning
    output:
    set sid, val(bname), "${sid}__*_cleaned.trk" optional true into bundle_for_eval
    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
    """
    scil_outlier_rejection.py ${bundle} "${sid}__${bname}_cleaned.trk" --alpha $outlier_alpha
    """
}

ref_bundles
    .transpose()
    // Get bundle name from ref bundles
    // (Remove everything before the first '__' and after the first '.' )
    .map{ it -> [it[0], it[1].name.replaceFirst(/.*__/, "").replaceFirst(/\..*/, ""), it[1]] }
    // Merge with bundle_for_eval using [sid,bname] as key
    .join(bundle_for_eval, by:[0,1], remainder: true)
    // Remove any candidate bundles not present in the gold standard
    .filter { it[2] != null }
    .branch {
        candidate_found: it[3] != null
            return it
        candidate_missing: it[3] == null
            return it
    }
    .set { branch_result }

branch_result.candidate_found
    .set { bundle_and_ref_for_eval }

branch_result.candidate_missing
    .set { bundle_missing_and_ref_for_eval }


process Compute_Measures {
    input:
    set sid, val(bname), file(ref_bundle), file(candidate_bundle) from bundle_and_ref_for_eval
    output:
    set sid, "${sid}__${bname}_individual_measures.json", "${sid}__${bname}_pairwise_measures.json" into bundle_measures_found
    script:
    """
    scil_evaluate_bundles_individual_measures.py ${candidate_bundle} ${sid}__${bname}_individual_measures.json
    evaluate_candidate_bundle.py --in ${candidate_bundle} --gs ${ref_bundle} --out ${sid}__${bname}_pairwise_measures.json
    """
}

process Compute_Measures_Missing {
    input:
    set sid, val(bname), file(ref_bundle), val(candidate_bundle) from bundle_missing_and_ref_for_eval
    output:
    set sid, "${sid}__${bname}_individual_measures.json", "${sid}__${bname}_pairwise_measures.json" into bundle_measures_missing
    script:
//  If candidate bundle is absent, return empty JSON file.
    """
    echo "{}" > ${sid}__${bname}_individual_measures.json
    echo "{}" > ${sid}__${bname}_pairwise_measures.json
    """
}

bundle_measures_found
    .mix(bundle_measures_missing)
    .collect {it -> [it[1], it[2]]}
    .set{all_bundles_measures}

process Aggregate_Measures {
    publishDir = "$params.output_dir/Aggregate_Measures"
    input:
    file(measures) from all_bundles_measures
    output:
    file "measures_aggregated.json"
    script:
    """
    aggregate_measures.py --in $measures --out measures_aggregated.json
    """
}
