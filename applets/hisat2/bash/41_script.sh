#!/bin/bash

main() {

    # Download the inputs to the worker's local storage
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

    # Extract the tarball containing the HISAT2 reference
    tar xf hisat2_index.tar.gz

    # Get the index basename to use when calling hisat2
    index_filename=$(dx describe --name "$hisat2_index_targz")
    index_basename=${index_filename%.tar.gz}/genome

    # Run hisat2
    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    # Upload the resulting file and associate it with the aligned_sam output
    uploaded_id=$(dx upload hisat2_output.sam --brief)
    dx-jobutil-add-output aligned_sam $uploaded_id
}
