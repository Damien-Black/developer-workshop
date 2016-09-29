#!/bin/bash
set -x -o pipefail

main() {

    # Download the inputs to the worker's local storage
    dx-download-all-inputs --parallel
    gunzip "$mate1_fastqgz_path"
    mate1_fastq_path="${mate1_fastqgz_path%.gz}"
    mate1_fastq_name="${mate1_fastqgz_name%.gz}"
    gunzip "$mate2_fastqgz_path"
    mate2_fastq_path="${mate2_fastqgz_path%.gz}"
    mate2_fastq_name="${mate2_fastqgz_name%.gz}"

    # Extract the tarball containing the HISAT2 reference
    mark-section "extracting reference tarball"
    tar xf "$hisat2_index_targz_path"

    # Get the index basename to use when calling hisat2
    index_directory=${hisat2_index_targz_name%.tar.gz}
    index1_filename=( $index_directory/*.1.ht2 )
    index_basename=${index1_filename%.1.ht2}

    # Run hisat2
    mark-section "running hisat2 command"
    hisat2 --dta -x $index_basename -1 "$mate1_fastq_path" -2 "$mate2_fastq_path" -S hisat2_output.sam

    mark-section "converting SAM to BAM"
    mkdir -p ./out/aligned_bam
    samtools view -b hisat2_output.sam > ./out/aligned_bam/"$mate1_fastqgz_prefix".hisat2.bam

    dx-upload-all-outputs

    mark-success
}
