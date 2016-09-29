import glob
import multiprocessing
import os
import subprocess
import traceback

import dxpy

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastqgz, mate2_fastqgz):

    # First, download all the input files to local storage
    inputs_dict = dxpy.download_all_inputs()
    
    gunzip_proc_1 = subprocess.Popen(["gunzip", inputs_dict["mate1_fastqgz_path"][0]])
    gunzip_proc_2 = subprocess.Popen(["gunzip", inputs_dict["mate2_fastqgz_path"][0]])

    gunzip_proc_1.wait()
    gunzip_proc_2.wait()

    if any([p.returncode != 0 for p in [gunzip_proc_1, gunzip_proc_2]]):
        raise dxpy.AppError("Gunzipping FASTQ files failed.")


    # Second, extract the index tarball
    try:
        subprocess.check_call(["tar", "xf", inputs_dict["hisat2_index_targz_path"][0]])
    except Exception:
        traceback.print_exc()
        raise dxpy.AppError("Error extracting reference tarball")

    index_directory = inputs_dict["hisat2_index_targz_name"][0][:-len(".tar.gz")]
    index_basename = glob.glob(os.path.join(index_directory, "*.1.ht2"))[0][:-len(".1.ht2")]


    # Prepare the hisat2 command and run it.
    hisat2_cmd_template = ("hisat2 --dta -p {nproc} -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq}")
    hisat2_cmd = hisat2_cmd_template.format(
        nproc=multiprocessing.cpu_count(),
        index_basename=index_basename,
        mate1_fastq=inputs_dict["mate1_fastqgz_path"][0][:-len(".gz")],
        mate2_fastq=inputs_dict["mate2_fastqgz_path"][0][:-len(".gz")])
    hisat2_proc = subprocess.Popen(hisat2_cmd, shell=True, stdout=subprocess.PIPE)

    samtools_proc = subprocess.Popen(
        ["samtools", "view", "-b", "-",
         "-o", str(inputs_dict["mate1_fastqgz_prefix"][0]) + ".bam"],
        stdin=hisat2_proc.stdout)
    hisat2_proc.stdout.close()
    samtools_proc.communicate()
    hisat2_proc.wait()

    if hisat2_proc.returncode != 0 or samtools_proc.returncode != 0:
        raise dxpy.AppError("hisat2 | samtools pipeline failed")

    # Upload the output SAM file.
    uploaded_dxfile = dxpy.upload_local_file(
        inputs_dict["mate1_fastqgz_prefix"][0] + ".bam")

    # Return the ID of the uploaded SAM file associated with the "aligned_sam"
    # field in the outputSpec in dxapp.json.
    return {"aligned_bam": dxpy.dxlink(uploaded_dxfile.get_id())}
