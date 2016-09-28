---
title: Improving File IO
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_improving_file_io.html
---

If unmanaged, file IO can be a significant source of overhead when running
jobs on DNAnexus. But, there are a number of ways to ensure the file transfers
are efficient and easy for users.

## 7.0 File Patterns

When you are selecting files as inputs through the web UI, you can limit the files that are
suggested to files that match a certain pattern. For example, for the FASTQ file inputs, you
can limit suggested files to those that end in ".fastq" or ".fq". For the reference tarball,
you can limit files to those ending with ".tar.gz". This is done with the "patterns" key in
the dxapp.json file:

```json
{
  "name": "hisat2",
  "title": "HISAT2",
  "summary": "Runs a simple HISAT2 command.",
  "version": "0.0.1",
  "inputSpec": [
    {
      "name": "hisat2_index_targz",
      "label": "HISAT2 Index Tarball",
      "class": "file",
      "help": "Tar.gz'd HISAT2 index for a reference genome, produced using hisat2-build.",
      "patterns": ["*.tar.gz"]
    },
    {
      "name": "mate1_fastq",
      "label": "Mate 1 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 1 reads.",
      "patterns": ["*.fastq", "*.fq"]
    },
    {
      "name": "mate2_fastq",
      "label": "Mate 2 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 2 reads.",
      "patterns": ["*.fastq", "*.fq"]
    }
  ],
  "outputSpec": [
    {
      "name": "aligned_sam",
      "label": "Aligned SAM",
      "class": "file",
      "help": "SAM file with alignments reported by HISAT2",
      "patterns": ["*.sam"]
    }
  ],
  "runSpec": {
    "file": "src/script.sh",
    "interpreter": "bash",
    "systemRequirements": {"*": {"instanceType": "mem1_ssd1_x8"}},
    "distribution": "Ubuntu",
    "release": "14.04",
    "timeoutPolicy": {"*": {"hours": 6, "minutes": 30}}
  }
}
```

## 7.1 Compressing Inputs and Outputs

Generally it is faster to transfer a compressed file and decompress it than it is
to transfer an uncompressed file. In the current iteration of the HISAT2 applet,
our inputs are FASTQ files, and our output is a SAM. We can speed up the applet
by changing those to gzipped FASTQ files and a BAM file.

To do this, we need to change the dxapp.json file so that the patterns and names
now indicate that the user should expect compressed files:

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old" data-toggle="tab">old</a></li>
    <li><a href="#new" data-toggle="tab">new</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old">

<pre class="highlight"><code><span class="p">{</span><span class="w">
  </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"hisat2"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"title"</span><span class="p">:</span><span class="w"> </span><span class="s2">"HISAT2"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"summary"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Runs a simple HISAT2 command."</span><span class="p">,</span><span class="w">
  </span><span class="nt">"version"</span><span class="p">:</span><span class="w"> </span><span class="s2">"0.0.1"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"inputSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"hisat2_index_targz"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"HISAT2 Index Tarball"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Tar.gz'd HISAT2 index for a reference genome, produced using hisat2-build."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.tar.gz"</span><span class="p">]</span><span class="w">
    </span><span class="p">},</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mate1_fastq"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Mate 1 FASTQ"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"FASTQ file containing mate 1 reads."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.fastq"</span><span class="p">,</span><span class="w"> </span><span class="s2">"*.fq"</span><span class="p">]</span><span class="w">
    </span><span class="p">},</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mate2_fastq"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Mate 2 FASTQ"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"FASTQ file containing mate 2 reads."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.fastq"</span><span class="p">,</span><span class="w"> </span><span class="s2">"*.fq"</span><span class="p">]</span><span class="w">
    </span><span class="p">}</span><span class="w">
  </span><span class="p">],</span><span class="w">
  </span><span class="nt">"outputSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"aligned_sam"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Aligned SAM"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"SAM file with alignments reported by HISAT2"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.sam"</span><span class="p">]</span><span class="w">
    </span><span class="p">}</span><span class="w">
  </span><span class="p">],</span><span class="w">
  </span><span class="nt">"runSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="w">
    </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.sh"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"bash"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"systemRequirements"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"*"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"instanceType"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mem1_ssd1_x8"</span><span class="p">}},</span><span class="w">
    </span><span class="nt">"distribution"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Ubuntu"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"release"</span><span class="p">:</span><span class="w"> </span><span class="s2">"14.04"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"timeoutPolicy"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"*"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"hours"</span><span class="p">:</span><span class="w"> </span><span class="mi">6</span><span class="p">,</span><span class="w"> </span><span class="nt">"minutes"</span><span class="p">:</span><span class="w"> </span><span class="mi">30</span><span class="p">}}</span><span class="w">
  </span><span class="p">}</span><span class="w">
</span><span class="p">}</span><span class="w">
</span></code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="new">

<pre class="highlight"><code><span class="p">{</span><span class="w">
  </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"hisat2"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"title"</span><span class="p">:</span><span class="w"> </span><span class="s2">"HISAT2"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"summary"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Runs a simple HISAT2 command."</span><span class="p">,</span><span class="w">
  </span><span class="nt">"version"</span><span class="p">:</span><span class="w"> </span><span class="s2">"0.0.1"</span><span class="p">,</span><span class="w">
  </span><span class="nt">"inputSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"hisat2_index_targz"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"HISAT2 Index Tarball"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Tar.gz'd HISAT2 index for a reference genome, produced using hisat2-build."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.tar.gz"</span><span class="p">]</span><span class="w">
    </span><span class="p">},</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mate1_fastqgz"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Mate 1 FASTQ (gzipped)"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Gzipped FASTQ file containing mate 1 reads."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.fastq.gz"</span><span class="p">,</span><span class="w"> </span><span class="s2">"*.fq.gz"</span><span class="p">]</span><span class="w">
    </span><span class="p">},</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mate2_fastqgz"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Mate 2 FASTQ (gzipped)"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Gzipped FASTQ file containing mate 2 reads."</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.fastq.gz"</span><span class="p">,</span><span class="w"> </span><span class="s2">"*.fq.gz"</span><span class="p">]</span><span class="w">
    </span><span class="p">}</span><span class="w">
  </span><span class="p">],</span><span class="w">
  </span><span class="nt">"outputSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="w">
    </span><span class="p">{</span><span class="w">
      </span><span class="nt">"name"</span><span class="p">:</span><span class="w"> </span><span class="s2">"aligned_bam"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"label"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Aligned BAM"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"class"</span><span class="p">:</span><span class="w"> </span><span class="s2">"file"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"help"</span><span class="p">:</span><span class="w"> </span><span class="s2">"BAM file with alignments reported by HISAT2"</span><span class="p">,</span><span class="w">
      </span><span class="nt">"patterns"</span><span class="p">:</span><span class="w"> </span><span class="p">[</span><span class="s2">"*.bam"</span><span class="p">]</span><span class="w">
    </span><span class="p">}</span><span class="w">
  </span><span class="p">],</span><span class="w">
  </span><span class="nt">"runSpec"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="w">
    </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.sh"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"bash"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"systemRequirements"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"*"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"instanceType"</span><span class="p">:</span><span class="w"> </span><span class="s2">"mem1_ssd1_x8"</span><span class="p">}},</span><span class="w">
    </span><span class="nt">"distribution"</span><span class="p">:</span><span class="w"> </span><span class="s2">"Ubuntu"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"release"</span><span class="p">:</span><span class="w"> </span><span class="s2">"14.04"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"timeoutPolicy"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"*"</span><span class="p">:</span><span class="w"> </span><span class="p">{</span><span class="nt">"hours"</span><span class="p">:</span><span class="w"> </span><span class="mi">6</span><span class="p">,</span><span class="w"> </span><span class="nt">"minutes"</span><span class="p">:</span><span class="w"> </span><span class="mi">30</span><span class="p">}}</span><span class="w">
  </span><span class="p">}</span><span class="w">
</span><span class="p">}</span><span class="w">
</span></code></pre>
  
  </div>
</div>

And we need to update our applet scripts so that they decompress the input files and compress the output:

```shell
#!/bin/bash
set -x -o pipefail

main() {

    # Download the inputs to the worker's local storage
    mark-section "downloading input files"
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastqgz" -o mate1.fastq.gz
    gunzip mate1.fastq.gz
    dx download "$mate2_fastqgz" -o mate2.fastq.gz
    gunzip mate2.fastq.gz

    # Extract the tarball containing the HISAT2 reference
    mark-section "extracting reference tarball"
    tar xf hisat2_index.tar.gz

    # Get the index basename to use when calling hisat2
    index_filename=$(dx describe --name "$hisat2_index_targz")
    index_directory=${index_filename%.tar.gz}
    index1_filename=( $index_directory/*.1.ht2 )
    index_basename=${index1_filename%.1.ht2}

    # Run hisat2
    mark-section "running hisat2 command"
    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    mark-section "converting SAM to BAM"
    samtools view -b hisat2_output.sam > hisat2_output.bam

    # Upload the resulting file and associate it with the aligned_bam output
    mark-section "uploading BAM results"
    uploaded_id=$(dx upload hisat2_output.bam --brief)
    dx-jobutil-add-output aligned_bam "$uploaded_id"
    mark-success
}
```

```python
import glob
import os
import subprocess
import traceback

import dxpy

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastqgz, mate2_fastqgz):

    # First, download all the input files to local storage
    dxpy.download_dxfile(hisat2_index_targz, "hisat2_index.tar.gz")
    dxpy.download_dxfile(mate1_fastqgz, "mate1.fastq.gz")
    subprocess.check_call(["gunzip", "mate1.fastq.gz"])
    dxpy.download_dxfile(mate2_fastqgz, "mate2.fastq.gz")
    subprocess.check_call(["gunzip", "mate2.fastq.gz"])


    # Second, extract the index tarball
    try:
        subprocess.check_call(["tar", "xf", "hisat2_index.tar.gz"])
    except Exception:
        traceback.print_exc()
        raise dxpy.AppError("Error extracting reference tarball")

    index_directory = dxpy.DXFile(hisat2_index_targz).name[:-len(".tar.gz")]
    index_basename = glob.glob(os.path.join(index_directory, "*.1.ht2"))[0][:-len(".1.ht2")]


    # Prepare the hisat2 command and run it.
    hisat2_cmd_template = ("hisat2 --dta -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq} -S {hisat2_output_sam}")
    hisat2_cmd = hisat2_cmd_template.format(
        index_basename=index_basename,
        mate1_fastq="mate1.fastq",
        mate2_fastq="mate2.fastq",
        hisat2_output_sam="hisat2_output.sam")
    subprocess.check_call(hisat2_cmd, shell=True)

    subprocess.check_call(["samtools", "view", "-b", "hisat2_output.sam", "-o", "hisat2_output.bam"])

    # Upload the output SAM file.
    uploaded_dxfile = dxpy.upload_local_file("hisat2_output.bam")

    # Return the ID of the uploaded SAM file associated with the "aligned_sam"
    # field in the outputSpec in dxapp.json.
    return {"aligned_bam": dxpy.dxlink(uploaded_dxfile.get_id())}
```

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#python" data-toggle="tab">python</a></li>
    <li><a href="#bash" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="python">
<pre class="highlight"><code><span class="kn">import</span> <span class="nn">glob</span>
<span class="kn">import</span> <span class="nn">os</span>
<span class="kn">import</span> <span class="nn">subprocess</span>
<span class="kn">import</span> <span class="nn">traceback</span>

<span class="kn">import</span> <span class="nn">dxpy</span>

<span class="nd">@dxpy.entry_point</span><span class="p">(</span><span class="s">"main"</span><span class="p">)</span>
<span class="k">def</span> <span class="nf">main</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">,</span> <span class="n">mate1_fastqgz</span><span class="p">,</span> <span class="n">mate2_fastqgz</span><span class="p">):</span>

    <span class="c"># First, download all the input files to local storage</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">,</span> <span class="s">"hisat2_index.tar.gz"</span><span class="p">)</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">mate1_fastqgz</span><span class="p">,</span> <span class="s">"mate1.fastq.gz"</span><span class="p">)</span>
    <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">([</span><span class="s">"gunzip"</span><span class="p">,</span> <span class="s">"mate1.fastq.gz"</span><span class="p">])</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">mate2_fastqgz</span><span class="p">,</span> <span class="s">"mate2.fastq.gz"</span><span class="p">)</span>
    <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">([</span><span class="s">"gunzip"</span><span class="p">,</span> <span class="s">"mate2.fastq.gz"</span><span class="p">])</span>


    <span class="c"># Second, extract the index tarball</span>
    <span class="k">try</span><span class="p">:</span>
        <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">([</span><span class="s">"tar"</span><span class="p">,</span> <span class="s">"xf"</span><span class="p">,</span> <span class="s">"hisat2_index.tar.gz"</span><span class="p">])</span>
    <span class="k">except</span> <span class="nb">Exception</span><span class="p">:</span>
        <span class="n">traceback</span><span class="o">.</span><span class="n">print_exc</span><span class="p">()</span>
        <span class="k">raise</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">AppError</span><span class="p">(</span><span class="s">"Error extracting reference tarball"</span><span class="p">)</span>

    <span class="n">index_directory</span> <span class="o">=</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">DXFile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">)</span><span class="o">.</span><span class="n">name</span><span class="p">[:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".tar.gz"</span><span class="p">)]</span>
    <span class="n">index_basename</span> <span class="o">=</span> <span class="n">glob</span><span class="o">.</span><span class="n">glob</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">index_directory</span><span class="p">,</span> <span class="s">"*.1.ht2"</span><span class="p">))[</span><span class="mi">0</span><span class="p">][:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".1.ht2"</span><span class="p">)]</span>


    <span class="c"># Prepare the hisat2 command and run it.</span>
    <span class="n">hisat2_cmd_template</span> <span class="o">=</span> <span class="p">(</span><span class="s">"hisat2 --dta -x {index_basename} -1 {mate1_fastq} "</span>
                           <span class="s">"-2 {mate2_fastq} -S {hisat2_output_sam}"</span><span class="p">)</span>
    <span class="n">hisat2_cmd</span> <span class="o">=</span> <span class="n">hisat2_cmd_template</span><span class="o">.</span><span class="n">format</span><span class="p">(</span>
        <span class="n">index_basename</span><span class="o">=</span><span class="n">index_basename</span><span class="p">,</span>
        <span class="n">mate1_fastq</span><span class="o">=</span><span class="s">"mate1.fastq"</span><span class="p">,</span>
        <span class="n">mate2_fastq</span><span class="o">=</span><span class="s">"mate2.fastq"</span><span class="p">,</span>
        <span class="n">hisat2_output_sam</span><span class="o">=</span><span class="s">"hisat2_output.sam"</span><span class="p">)</span>
    <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">(</span><span class="n">hisat2_cmd</span><span class="p">,</span> <span class="n">shell</span><span class="o">=</span><span class="bp">True</span><span class="p">)</span>

    <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">([</span><span class="s">"samtools"</span><span class="p">,</span> <span class="s">"view"</span><span class="p">,</span> <span class="s">"-b"</span><span class="p">,</span> <span class="s">"hisat2_output.sam"</span><span class="p">,</span> <span class="s">"-o"</span><span class="p">,</span> <span class="s">"hisat2_output.bam"</span><span class="p">])</span>

    <span class="c"># Upload the output SAM file.</span>
    <span class="n">uploaded_dxfile</span> <span class="o">=</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">upload_local_file</span><span class="p">(</span><span class="s">"hisat2_output.bam"</span><span class="p">)</span>

    <span class="c"># Return the ID of the uploaded SAM file associated with the "aligned_sam"</span>
    <span class="c"># field in the outputSpec in dxapp.json.</span>
    <span class="k">return</span> <span class="p">{</span><span class="s">"aligned_bam"</span><span class="p">:</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">dxlink</span><span class="p">(</span><span class="n">uploaded_dxfile</span><span class="o">.</span><span class="n">get_id</span><span class="p">())}</span>
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="bash">

<pre class="highlight"><code><span class="c">#!/bin/bash</span>
<span class="nb">set</span> -x -o pipefail

main<span class="o">()</span> <span class="o">{</span>

    <span class="c"># Download the inputs to the worker's local storage</span>
    mark-section <span class="s2">"downloading input files"</span>
    dx download <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span> -o hisat2_index.tar.gz
    dx download <span class="s2">"</span><span class="nv">$mate1_fastqgz</span><span class="s2">"</span> -o mate1.fastq.gz
    gunzip mate1.fastq.gz
    dx download <span class="s2">"</span><span class="nv">$mate2_fastqgz</span><span class="s2">"</span> -o mate2.fastq.gz
    gunzip mate2.fastq.gz

    <span class="c"># Extract the tarball containing the HISAT2 reference</span>
    mark-section <span class="s2">"extracting reference tarball"</span>
    tar xf hisat2_index.tar.gz

    <span class="c"># Get the index basename to use when calling hisat2</span>
    <span class="nv">index_filename</span><span class="o">=</span><span class="k">$(</span>dx describe --name <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span><span class="k">)</span>
    <span class="nv">index_directory</span><span class="o">=</span><span class="k">${</span><span class="nv">index_filename</span><span class="p">%.tar.gz</span><span class="k">}</span>
    <span class="nv">index1_filename</span><span class="o">=(</span> <span class="nv">$index_directory</span>/<span class="k">*</span>.1.ht2 <span class="o">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index1_filename</span><span class="p">%.1.ht2</span><span class="k">}</span>

    <span class="c"># Run hisat2</span>
    mark-section <span class="s2">"running hisat2 command"</span>
    hisat2 --dta -x <span class="nv">$index_basename</span> -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    mark-section <span class="s2">"converting SAM to BAM"</span>
    samtools view -b hisat2_output.sam &gt; hisat2_output.bam

    <span class="c"># Upload the resulting file and associate it with the aligned_bam output</span>
    mark-section <span class="s2">"uploading BAM results"</span>
    <span class="nv">uploaded_id</span><span class="o">=</span><span class="k">$(</span>dx upload hisat2_output.bam --brief<span class="k">)</span>
    dx-jobutil-add-output aligned_bam <span class="s2">"</span><span class="nv">$uploaded_id</span><span class="s2">"</span>
    mark-success
<span class="o">}</span>
</code></pre>
  
  </div>
</div>

Note that we are now using `samtools` in out applet scripts, so we need to make sure
that is available to the worker. We can do that by either copying a samtools executable
into the applet's resources directory, or we can add it to the dxapp.json:


```json
  ...
  "runSpec": {
    ...
    "execDepends": [{"name": "samtools"}]
  }
  ...
```

This will make the worker run

```shell
apt-get install samtools
```
before the job starts.

## 7.2 Simplifying File IO Commands

### Bash

So far, we have been handling transfers to and from the platform with an individual command
for each file. We can make this easier by using a couple of commands that handle IO for
all input and output files at once. The first command is `dx-download-all-inputs`, and it
can replace the individual `dx download` commands we have been using:

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old2" data-toggle="tab">old</a></li>
    <li><a href="#new2" data-toggle="tab">new</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old2">

<pre class="highlight"><code>    dx download <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span> -o hisat2_index.tar.gz
    dx download <span class="s2">"</span><span class="nv">$mate1_fastqgz</span><span class="s2">"</span> -o mate1.fastq.gz
    dx download <span class="s2">"</span><span class="nv">$mate2_fastqgz</span><span class="s2">"</span> -o mate2.fastq.gz
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="new2">

<pre class="highlight"><code>   dx-download-all-inputs
</code></pre>
  
  </div>
</div>

Now, we no longer have hard-coded file names for the inputs that we have downloaded, but we still
need to refer to them later in the script. We can do this using special variables that are
available in the bash script. For example, supposed the name of the file we were supplying for the
`mate1_fastqgz` file was `SRR3469431_1.fastq.gz`. Then, the three variables we can use in the bash
script would have the following values:

```shell
echo "$mate1_fastqgz_path"   # /home/dnanexus/in/mate1_fastqgz/SRR3469431_1.fastq.gz 
echo "$mate1_fastqgz_name"   # SRR3469431_1.fastq.gz 
echo "$mate1_fastqgz_prefix" # SRR3469431_1
```

So, we can just refer to the downloaded file using these variables. If the input were an
array of files rather than a single file, the correponding bash variables would be arrays
rather than strings.

`dx-upload-all-outputs` handles both the uploading of files and the association of the
uploaded files with the appropriate fields in the dxapp.json's outputSpec. To use
`dx-upload-all-outputs`, move output files to directories `/home/dnanexus/<output_field>/`.
So for the `aligned_bam` output, we should put our BAM file in `/home/dnanexuss/aligned_bam/`
before calling `dx-upload-all-outputs`.

Putting these together, we get an updated script:

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old3" data-toggle="tab">old</a></li>
    <li><a href="#new3" data-toggle="tab">new</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old3">

<pre class="highlight"><code><span class="c">#!/bin/bash</span>
<span class="nb">set</span> -x -o pipefail

main<span class="o">()</span> <span class="o">{</span>

    <span class="c"># Download the inputs to the worker's local storage</span>
    mark-section <span class="s2">"downloading input files"</span>
    dx download <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span> -o hisat2_index.tar.gz
    dx download <span class="s2">"</span><span class="nv">$mate1_fastqgz</span><span class="s2">"</span> -o mate1.fastq.gz
    gunzip mate1.fastq.gz
    dx download <span class="s2">"</span><span class="nv">$mate2_fastqgz</span><span class="s2">"</span> -o mate2.fastq.gz
    gunzip mate2.fastq.gz

    <span class="c"># Extract the tarball containing the HISAT2 reference</span>
    mark-section <span class="s2">"extracting reference tarball"</span>
    tar xf hisat2_index.tar.gz

    <span class="c"># Get the index basename to use when calling hisat2</span>
    <span class="nv">index_filename</span><span class="o">=</span><span class="k">$(</span>dx describe --name <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span><span class="k">)</span>
    <span class="nv">index_directory</span><span class="o">=</span><span class="k">${</span><span class="nv">index_filename</span><span class="p">%.tar.gz</span><span class="k">}</span>
    <span class="nv">index1_filename</span><span class="o">=(</span> <span class="nv">$index_directory</span>/<span class="k">*</span>.1.ht2 <span class="o">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index1_filename</span><span class="p">%.1.ht2</span><span class="k">}</span>

    <span class="c"># Run hisat2</span>
    mark-section <span class="s2">"running hisat2 command"</span>
    hisat2 --dta -x <span class="nv">$index_basename</span> -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    mark-section <span class="s2">"converting SAM to BAM"</span>
    samtools view -b hisat2_output.sam &gt; hisat2_output.bam

    <span class="c"># Upload the resulting file and associate it with the aligned_bam output</span>
    mark-section <span class="s2">"uploading BAM results"</span>
    <span class="nv">uploaded_id</span><span class="o">=</span><span class="k">$(</span>dx upload hisat2_output.bam --brief<span class="k">)</span>
    dx-jobutil-add-output aligned_bam <span class="s2">"</span><span class="nv">$uploaded_id</span><span class="s2">"</span>
    mark-success
<span class="o">}</span>
</code></pre>


  </div>
  <div role="tabpanel" class="tab-pane" id="new3">

<pre class="highlight"><code><span class="c">#!/bin/bash</span>
<span class="nb">set</span> -x -o pipefail

main<span class="o">()</span> <span class="o">{</span>

    <span class="c"># Download the inputs to the worker's local storage</span>
    dx-download-all-inputs --parallel
    gunzip <span class="s2">"</span><span class="nv">$mate1_fastqgz_path</span><span class="s2">"</span>
    <span class="nv">mate1_fastq_path</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">mate1_fastqgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>
    <span class="nv">mate1_fastq_name</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">mate1_fastqgz_name</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>
    gunzip <span class="s2">"</span><span class="nv">$mate2_fastqgz_path</span><span class="s2">"</span>
    <span class="nv">mate2_fastq_path</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">mate2_fastqgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>
    <span class="nv">mate2_fastq_name</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">mate2_fastqgz_name</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>

    <span class="c"># Extract the tarball containing the HISAT2 reference</span>
    mark-section <span class="s2">"extracting reference tarball"</span>
    tar xf <span class="s2">"</span><span class="nv">$hisat2_index_targz_path</span><span class="s2">"</span>

    <span class="c"># Get the index basename to use when calling hisat2</span>
    <span class="nv">index_directory</span><span class="o">=</span><span class="k">${</span><span class="nv">hisat2_index_targz_name</span><span class="p">%.tar.gz</span><span class="k">}</span>
    <span class="nv">index1_filename</span><span class="o">=(</span> <span class="nv">$index_directory</span>/<span class="k">*</span>.1.ht2 <span class="o">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index1_filename</span><span class="p">%.1.ht2</span><span class="k">}</span>

    <span class="c"># Run hisat2</span>
    mark-section <span class="s2">"running hisat2 command"</span>
    hisat2 --dta -x <span class="nv">$index_basename</span> -1 <span class="s2">"</span><span class="nv">$mate1_fastq_path</span><span class="s2">"</span> -2 <span class="s2">"</span><span class="nv">$mate2_fastq_path</span><span class="s2">"</span> -S hisat2_output.sam

    mark-section <span class="s2">"converting SAM to BAM"</span>
    mkdir -p ./out/aligned_bam
    samtools view -b hisat2_output.sam &gt; ./out/aligned_bam/<span class="s2">"</span><span class="nv">$mate1_fastqgz_prefix</span><span class="s2">"</span>.hisat2.bam

    dx-upload-all-outputs

    mark-success
<span class="o">}</span>
</code></pre>
  
  </div>
</div>

### Python

There is no equivalent for `dx-upload-all-outputs` in python, but there is a function
`dxpy.dx_download_all_inputs()`. This will return a dictionary with keys that are the
same as the variables used in bash.


## 7.3 Piping

So far, we have been completely downloading each compressed file and then decompressing it.
We can speed things up by decompressing as we download, so that the file is streamed to a
`gunzip` process as it is downloaded. We do this by using pipes and `dx cat`:

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old4" data-toggle="tab">old</a></li>
    <li><a href="#new4" data-toggle="tab">new</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old4">

<pre class="highlight"><code>    dx download <span class="s2">"</span><span class="nv">$mate1_fastqgz</span><span class="s2">"</span> -o mate1.fastq.gz
    gunzip mate1.fastq.gz
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="new4">

<pre class="highlight"><code>   dx cat <span class="s2">"</span><span class="nv">$mate1_fastqgz</span><span class="s2">"</span> | gunzip &gt; mate1_fastqgz
</code></pre>
  
  </div>
</div>

We also would like to keep using the file variables, which takes a little special handling:

```shell
    mate1_fastq_path="${mate1_fastqgz_path%.gz}"
    mkdir -p "$(dirname "$mate1_fastq_path")"
    dx cat "$mate1_fastqgz" | gunzip > "$mate1_fastq_path"
```

In python this can be more complicated. One option is to just use `subprocess` with `shell=True`, and
simply run commands the way you would in bash. Alternatively, you can use `subprocess.PIPE`:

```python
    # Prepare the hisat2 command and run it.
    hisat2_cmd_template = ("hisat2 --dta -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq}")
    hisat2_cmd = hisat2_cmd_template.format(
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
```

In this code, the `subprocess.Popen` call for the hisat2 command sends its output
to the `samtools_proc` process. This can be tricky to get right, so it is helpful to
start from a working template.
