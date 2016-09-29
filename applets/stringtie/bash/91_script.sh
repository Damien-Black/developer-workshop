#!/bin/bash
set -x -o pipefail

stringtie_on_one_bam() {
    dx-download-all-inputs
    gunzip "$gtfgz_path"
    samtools sort -@ "$(nproc)" "$bam_path" -o sorted.bam
    mv sorted.bam "$bam_path"
    mkdir -p ./out/assembled_transcripts_gtf
    stringtie "$bam_path" -o ./out/assembled_transcripts_gtf/"$bam_prefix".assembled_transcripts.gtf \
        -p "$(nproc)" -G "${gtfgz_path%.gz}"

    dx-upload-all-outputs
}

stringtie_merge() {

    dx-download-all-inputs
    
    gunzip "$guide_gtfgz_path"
    mkdir -p ./out/merged_assembled_transcripts_gtf
    stringtie -G "${guide_gtfgz_path%.gz}" --merge "${assembled_transcripts_gtfs_path[@]}" \
        -o ./out/merged_assembled_transcripts_gtf/merged.assembled_transcripts.gtf -p "$(nproc)"
    
    dx-upload-all-outputs
}

stringtie_prepare_ballgown() {
    
    dx-download-all-inputs
    samtools sort -@ "$(nproc)" "$bam_path" -o sorted.bam
    mv sorted.bam "$bam_path"

    stringtie -eB -G "$merged_assembled_transcripts_gtf_path" "$bam_path" \
        -p "$(nproc)" -o /home/dnanexus/stringtie_output/output.gtf

    mkdir -p ./out/ballgown_ctabs_targz
    cd stringtie_output && tar czvf ./out/ballgown_ctabs_targz/"$bam_prefix".ballgown_ctabs.tar.gz *.ctab && cd ..
    dx-upload-all-outputs
}

main() {
    
    for bam in "${bams[@]}"; do
        one_bam_jobs+=($(dx-jobutil-new-job -ibam="$bam" -igtfgz="$reference_gtfgz" stringtie_on_one_bam))
    done
    
    for one_bam_job in "${one_bam_jobs[@]}"; do 
        merge_bam_args+=(-iassembled_transcripts_gtfs="$one_bam_job":assembled_transcripts_gtf)
    done
    merge_job=$(dx-jobutil-new-job -iguide_gtfgz="$reference_gtfgz" "${merge_bam_args[@]}" stringtie_merge)

    for bam in "${bams[@]}"; do
        ballgown_jobs+=($(dx-jobutil-new-job -ibam="$bam" \
            -imerged_assembled_transcripts_gtf="$merge_job":merged_assembled_transcripts_gtf \
            stringtie_prepare_ballgown))
    done
    
    for one_bam_job in "${one_bam_jobs[@]}"; do 
        dx-jobutil-add-output assembled_transcripts_gtfs \
            "$one_bam_job":assembled_transcripts_gtf --class=array:jobref
    done

    dx-jobutil-add-output merged_assembled_transcripts_gtf "$merge_job":merged_assembled_transcripts_gtf \
        --class=jobref

    for ballgown_job in "${ballgown_jobs[@]}"; do
        dx-jobutil-add-output ballgown_ctabs_targzs "$ballgown_job":ballgown_ctabs_targz --class=array:jobref
    done
}
