#!/bin/bash
set -x -o pipefail

main() {

    # Download the inputs to the worker's local storage
    mark-section "downloading inputs"

    mate1_fastq_path="${mate1_fastqgz_path%.gz}"
    mkdir -p "$(dirname "$mate1_fastq_path")"
    dx cat "$mate1_fastqgz" | gunzip > "$mate1_fastq_path"

    mate2_fastq_path="${mate2_fastqgz_path%.gz}"
    mkdir -p "$(dirname "$mate2_fastq_path")"
    dx cat "$mate2_fastqgz" | gunzip > "$mate2_fastq_path"

    index_directory=${hisat2_index_targz_name%.tar.gz}
    dx cat "$hisat2_index_targz" | tar xzvf - -C ~

    # Get the index basename to use when calling hisat2
    index1_filename=( $index_directory/*.1.ht2 )
    index_basename=${index1_filename%.1.ht2}

    # Run hisat2
    mark-section "running hisat2 command"
    hisat2 --dta -p "$(nproc)" -x $index_basename -1 "$mate1_fastq_path" -2 "$mate2_fastq_path" -S hisat2_output.sam

    mark-section "converting SAM to BAM"
    mkdir -p ./out/aligned_bam
    samtools view -b hisat2_output.sam > ./out/aligned_bam/"$mate1_fastqgz_prefix".hisat2.bam

    dx-upload-all-outputs

    mark-success
}
