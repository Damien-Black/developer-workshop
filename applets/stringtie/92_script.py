import multiprocessing
import shutil
import subprocess
import dxpy

@dxpy.entry_point("stringtie_on_one_bam")
def stringtie_on_one_bam(bam, gtfgz):
    bam_name = dxpy.DXFile(bam).name
    dxpy.download_dxfile(bam, bam_name)

    gtfgz_name = dxpy.DXFile(gtfgz).name
    dxpy.download_dxfile(gtfgz, gtfgz_name)

    subprocess.check_call(["gunzip", gtfgz_name])
    subprocess.check_call(["samtools", "sort", "-@", multiprocessing.cpu_count(),
                           bam_name, "-o", "sorted.bam"])
    shutil.move("sorted.bam", bam_name)
    
    output_gtf = bam_name[:-len(".bam")] + ".assembled_transcripts.gtf"
    stringtie_cmd = "stringtie {bam} -o {output_gtf} -p {nproc} -G {gtf}".format(
        bam=bam_name,
        output_gtf=output_gtf,
        nproc=multiprocessing.cpu_count(),
        gtf=gtfgz_name[:-len(".gz")])
    subprocess.check_call(stringtie_cmd, shell=True)

    uploaded_id = dxpy.upload_local_file(output_gtf).get_id()

    return {"assembled_transcripts_gtf": dxpy.dxlink(uploaded_id)} 

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
