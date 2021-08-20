#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_config":"$params.atlas_config",
                "atlas_directory":"$params.atlas_directory",
                "atlas_centroids":"$params.atlas_centroids",
                "multi_parameters":"$params.multi_parameters",
                "minimal_vote_ratio":"$params.minimal_vote_ratio",
                "wb_clustering_thr":"$params.wb_clustering_thr",
                "seeds":"$params.seeds",
                "outlier_alpha":"$params.outlier_alpha",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL TractoEval pipeline"
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
log.info "Atlas Config: $params.atlas_config"
log.info "Atlas Anat: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info "Atlas Centroids: $params.atlas_centroids"
log.info ""
log.info "[Recobundles options]"
log.info "Multi-Parameters Executions: $params.multi_parameters"
log.info "Minimal Vote Percentage: $params.minimal_vote_ratio"
log.info "Whole Brain Clustering Threshold: $params.wb_clustering_thr"
log.info "Random Seeds: $params.seeds"
log.info "Outlier Removal Alpha: $params.outlier_alpha"
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
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_centroids;anat_for_reference_bundles}

ref_bundles = Channel
    .fromFilePairs("$ref/**/{*.trk,}",
                    size:-1,
                    maxDepth:1) {it.parent.name}

if (!(params.atlas_anat) || !(params.atlas_config) || !(params.atlas_directory)) {
    error "You must specify all 3 atlas related input. --atlas_anat, " +
    "--atlas_config and --atlas_directory all are mandatory."
}

Channel.fromPath("$params.atlas_anat")
    .into{atlas_anat;atlas_anat_for_average}
atlas_config = Channel.fromPath("$params.atlas_config")
atlas_directory = Channel.fromPath("$params.atlas_directory")

if (params.atlas_centroids) {
    atlas_centroids = Channel.fromPath("$params.atlas_centroids/*_centroid.trk")
}
else {
    atlas_centroids = Channel.empty()
}

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
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition, transformation_for_centroids
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_average
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__native_anat.nii.gz"
    script:
//     """
//     export ANTS_RANDOM_SEED=1234
//     antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
//     cp ${native_anat} ${sid}__native_anat.nii.gz
//     """
    """
    touch ${sid}__output0GenericAffine.mat
    touch "${sid}__outputWarped.nii.gz"
    touch "${sid}__native_anat.nii.gz"
    """
}


anat_for_reference_centroids
    .join(transformation_for_centroids, by: 0)
    .set{anat_and_transformation}
process Transform_Centroids {
    input:
    set sid, file(anat), file(transfo) from anat_and_transformation
    each file(centroid) from atlas_centroids
    output:
    file "${sid}__${centroid.baseName}.trk"
    script:
//     """
//     scil_apply_transform_to_tractogram.py ${centroid} ${anat} ${transfo} ${sid}__${centroid.baseName}.trk --inverse --cut_invalid
//     """
        """
        touch "${sid}__${centroid.baseName}.trk"
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
    set sid, "*.trk" into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"
    script:
//     """
//     if [ `echo $tractograms | wc -w` -gt 1 ]; then
//         scil_streamlines_math.py lazy_concatenate $tractograms tracking_concat.trk
//     else
//         mv $tractograms tracking_concat.trk
//     fi
//     scil_remove_invalid_streamlines.py tracking_concat.trk tractogram_ic.trk --reference ${reference} --remove_single_point --remove_overlapping_points
//     mkdir tmp/
//     scil_recognize_multi_bundles.py tractogram_ic.trk ${config} ${directory}/*/ ${transfo} --inverse --out_dir tmp/ \
//         --log_level DEBUG --multi_parameters $params.multi_parameters --minimal_vote_ratio $params.minimal_vote_ratio \
//         --tractogram_clustering_thr $params.wb_clustering_thr --seeds $params.seeds --processes $params.rbx_processes
//     rm tractogram_ic.trk tracking_concat.trk
//     mv tmp/* ./
//     """
    """
    touch "results.json"
    touch "logfile.txt"
    touch AF_L.trk AF_R.trk CC_Oc.trk
    """
}

bundles_for_cleaning
    .transpose()
    .set{all_bundles_for_cleaning}
process Clean_Bundles {
    input:
    set sid, file(bundle) from all_bundles_for_cleaning
    output:
    set sid, val(bname), "${sid}__*_cleaned.trk" into bundle_for_eval
    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
//     """
//     scil_outlier_rejection.py ${bundle} "${sid}__${bname}_cleaned.trk" --alpha $params.outlier_alpha
//     """
    """
    touch "${sid}__${bname}_cleaned.trk"
    """
}

ref_bundles
    .view()
    .transpose()
    .view()
    // Get bundle name from ref bundles
    // (Remove everything before the first '__' and after the first '.' )
    .map{ it -> [it[0], it[1].name.replaceFirst(/.*__/, "").replaceFirst(/\..*/, ""), it[1]] }
    // Merge with bundle_for_eval using [sid,bname] as key
    .join(bundle_for_eval, by:[0,1])
    .view()
    .set{bundle_and_ref_for_eval}

process Compute_Measures {
    input:
    set sid, val(bname), file(ref_bundle), file(candidate_bundle) from bundle_and_ref_for_eval

    output:
    file "${sid}__${bname}__individual_measures.json"
    file "${sid}__${bname}__pairwise_measures.json"
    script:
    """
    scil_evaluate_bundles_individual_measures.py ${candidate_bundle} ${sid}__${bname}__individual_measures.json
    compute_pairwise_measures.py ${candidate_bundle} ${ref_bundle} ${sid}__${bname}__pairwise_measures.json
    """
}

// process Aggregate_Metrics {
//     input:
//
//     output:
//
//     script:
// }
