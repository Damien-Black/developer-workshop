import multiprocessing
import os
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
    subprocess.check_call(["samtools", "sort", "-@", str(multiprocessing.cpu_count()),
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

@dxpy.entry_point("stringtie_merge")
def stringtie_merge(guide_gtf, assembled_transcripts_gtfs):

    inputs_dict = dxpy.download_all_inputs()
    subprocess.check_call(["gunzip", inputs_dict["guide_gtf_path"][0]])
    stringtie_cmd = ("stringtie -G {guide_gtf} --merge {assembled_transcripts_gtfs} "
                     "-o {merged_output} -p {nproc}").format(
                             guide_gtf=inputs_dict["guide_gtf_path"][0][:-len(".gz")],
                             assembled_transcripts_gtfs=' '.join(inputs_dict["assembled_transcripts_gtfs_path"]),
                             merged_output="merged_assembled_transcripts.gtf",
                             nproc=multiprocessing.cpu_count())
    subprocess.check_call(stringtie_cmd, shell=True)

    uploaded_id = dxpy.upload_local_file("merged_assembled_transcripts.gtf").get_id()

    return {"merged_assembled_transcripts_gtf": dxpy.dxlink(uploaded_id)}

@dxpy.entry_point("stringtie_prepare_ballgown")
def stringtie_prepare_ballgown(bam, merged_assembled_transcripts_gtf):

    inputs_dict = dxpy.download_all_inputs()
    subprocess.check_call(["samtools", "sort", "-@", str(multiprocessing.cpu_count()),
                           inputs_dict["bam_path"][0], "-o", "sorted.bam"])
    shutil.move("sorted.bam", inputs_dict["bam_path"][0])
    os.makedirs("stringtie_output")
    
    stringtie_cmd = "stringtie -eB -G {merged_transcripts} {bam} -p {nproc} -o {output}".format(
        merged_transcripts=inputs_dict["merged_assembled_transcripts_gtf_path"][0],
        bam=inputs_dict["bam_path"][0],
        nproc=multiprocessing.cpu_count(),
        output=os.path.join("stringtie_output", inputs_dict["bam_prefix"][0] + '.gtf'))

    subprocess.check_call(stringtie_cmd, shell=True)

    subprocess.check_call(('cd stringtie_output && '
                           'tar czvf /home/dnanexus/' + inputs_dict["bam_prefix"][0] + '.ballgown_ctabs.tar.gz *.ctab && '
                           'cd ..'), shell=True)

    uploaded_id = dxpy.upload_local_file(inputs_dict["bam_prefix"][0] + '.ballgown_ctabs.tar.gz')

    return {"ballgown_ctabs_targz": dxpy.dxlink(uploaded_id)}

@dxpy.entry_point("main")
def main(bams, reference_gtfgz):
    
    one_bam_jobs = []
    for bam in bams:
        one_bam_job = dxpy.new_dxjob(
            fn_input={"bam": bam, "gtfgz": reference_gtfgz},
            fn_name="stringtie_on_one_bam")
        one_bam_jobs.append(one_bam_job)
            
    merge_job = dxpy.new_dxjob(
        fn_input={"guide_gtf": reference_gtfgz,
                    "assembled_transcripts_gtfs": [j.get_output_ref("assembled_transcripts_gtf") for j in one_bam_jobs]},
        fn_name="stringtie_merge")
    
    ballgown_jobs = []
    for bam in bams:
        ballgown_job = dxpy.new_dxjob(
            fn_input={"bam": bam, "merged_assembled_transcripts_gtf": merge_job.get_output_ref("merged_assembled_transcripts_gtf")},
            fn_name="stringtie_prepare_ballgown")
        ballgown_jobs.append(ballgown_job)

    return {"assembled_transcripts_gtfs": [dxpy.dxlink(j.get_output_ref("assembled_transcripts_gtf")) for j in one_bam_jobs],
            "merged_assembled_transcripts_gtf": dxpy.dxlink(merge_job.get_output_ref("merged_assembled_transcripts_gtf")),
            "ballgown_ctabs_targzs": [dxpy.dxlink(j.get_output_ref("ballgown_ctabs_targz")) for j in ballgown_jobs]}
