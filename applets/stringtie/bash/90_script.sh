#!/bin/bash
set -x -o pipefail

dx-download-all-inputs

gunzip "$reference_gtfgz_path"

mkdir -p ./out/assembled_transcripts_gtfs

for idx in "${!bams_path[@]}"; do
    bam_path="${bams_path[$idx]}"
    bam_prefix="${bams_prefix[$idx]}"

    samtools sort -@ "$(nproc)" "${bams_path[$idx]}" -o sorted.bam
    mv sorted.bam "${bams_path[$idx]}"

    stringtie "$bam_path" -o ./out/assembled_transcripts_gtfs/"$bam_prefix".assembled_transcripts.gtf \
        -p "$(nproc)" -G "${reference_gtfgz_path%.gz}"
    assembled_transcripts_gtfs+=(./out/assembled_transcripts_gtfs/"$bam_prefix".assembled_transcripts.gtf)
done

mkdir -p ./out/merged_assembled_transcripts_gtf
stringtie -G "${reference_gtfgz_path%.gz}" --merge "${assembled_transcripts_gtfs[@]}" \
    -o ./out/merged_assembled_transcripts_gtf/merged_assembled_transcripts.gtf -p "$(nproc)"

mkdir -p ./out/ballgown_ctabs_targzs
for idx in "${!bams_path[@]}"; do

    bam_path="${bams_path[$idx]}"
    bam_prefix="${bams_prefix[$idx]}"

    stringtie -eB -G ./out/merged_assembled_transcripts_gtf/merged_assembled_transcripts.gtf \
        "$bam_path" -p "$(nproc)" -o /home/dnanexus/"$bam_prefix"_stringtie_output/"$bam_prefix".gtf

    cd /home/dnanexus/"$bam_prefix"_stringtie_output/ && \
        tar czvf /home/dnanexus/out/ballgown_ctabs_targzs/"$bam_prefix".ballgown_ctabs.tar.gz *.ctab && \
        cd ..

done 

dx-upload-all-outputs
