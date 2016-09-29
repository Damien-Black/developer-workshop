---
title: Distributed Applets
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_distributed_applets.html
---

## 9.0 StringTie on a Single Worker

Now we are going to work on the StringTie applet to illustrate how to distribute work
over mutiple workers within one applet. 

Recall that the input to the StringTie applet will be a list of BAM files and a reference
GTF file. We will need to run three different commands:

1. For each BAM file
```shell
$ stringtie <bam> -G <reference_gtf> -o <assembled_transcripts>
```
2. Once for all assembled transcripts
```shell
stringtie -G <reference_gtf> --merge <assembled_transcripts_1> <assembled_transcripts_2> ... -o <merged_assembled_transcripts>
```
3. For each BAM file, create the `ctab` files for Ballgown
```shell
stringtie -eB -G <merged_assembled_transcripts> <bam> -o <output_gtf>
```

Using the practices we discussed when creating the HISAT2 applet, we can prepare the code for
an applet that runs all of these commands on a single DNAnexus worker.

```shell
#!/bin/bash
set -x -o pipefail

dx-download-all-inputs

gunzip "$reference_gtfgz_path"

mkdir -p ./out/assembled_transcripts_gtfs

for idx in "${!bams_path[@]}"; do
    bam_path="${bams_path[$idx]}"
    bam_prefix="${bams_prefix[$idx]}"

    stringtie "$bam_path" -o ./out/assembled_transcripts_gtfs/"$bam_prefix".assembled_transcripts.gtf \
        -p "$(nproc)" -G "${reference_gtfgz_path%.gz}"
    assembled_transcripts_gtfs+=(./out/assembled_transcripts_gtfs/"bam_prefix".assembled_transcripts.gtf)
done

mkdir -p ./out/merged_assembled_transcripts_gtf
stringtie -G "${reference_gtfgz_path%.gz}" --merge "${assembled_transcripts_gtfs[@]}" \
    -o ./out/merged_assembled_transcripts_gtf/merged.assembled_transcripts.gtf -p "$(nproc)"

for idx in "${!bams_path[@]}"; do

    bam_path="${bams_path[$idx]}"
    bam_prefix="${bams_prefix[$idx]}"

    ctab_output_path=/home/dnanexus/stringtie_output/"$bam_prefix".ctab
    ctab_output_paths+=("$ctab_output_path")
    stringtie -eB -G ./out/merged_assembled_transcripts_gtf/merged_assembled_transcripts.gtf \
        "$bam_path" -p "$(nproc)" -o /home/dnanexus/stringtie_output/"$bam_prefix".gtf
done 

mkdir -p ./out/ballgown_ctabs
mv stringtie_output/*.ctab ./out/ballgown_ctabs

dx-upload-all-outputs
```

## 9.1 Distributing in Bash

While this will work, we may want to take advantage of the scale available on the DNAnexus
platfrom and distribute the tasks over multiple workers. Rather than running each step on
the same worker, we can run each combination of BAM and command on a different worker. To
do this, we need to break the parts of the applet that we want to run on different workers
into separate functions outside of the `main` function. Since each function will be run on
a separate worker, each function will also need to handle file IO.


<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#single" data-toggle="tab">single</a></li>
    <li><a href="#parallel" data-toggle="tab">parallel</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="single">

<pre class="highlight"><code><span class="k">for </span>idx <span class="k">in</span> <span class="s2">"</span><span class="k">${</span><span class="p">!bams_path[@]</span><span class="k">}</span><span class="s2">"</span>; <span class="k">do
    </span><span class="nv">bam_path</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">bams_path</span><span class="p">[</span><span class="nv">$idx</span><span class="p">]</span><span class="k">}</span><span class="s2">"</span>
    <span class="nv">bam_prefix</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">bams_prefix</span><span class="p">[</span><span class="nv">$idx</span><span class="p">]</span><span class="k">}</span><span class="s2">"</span>

    stringtie <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> -o ./out/assembled_transcripts_gtfs/<span class="s2">"</span><span class="nv">$bam_prefix</span><span class="s2">"</span>.assembled_transcripts.gtf <span class="se">\</span>
        -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> -G <span class="s2">"</span><span class="k">${</span><span class="nv">reference_gtfgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>
    assembled_transcripts_gtfs+<span class="o">=(</span>./out/assembled_transcripts_gtfs/<span class="s2">"bam_prefix"</span>.assembled_transcripts.gtf<span class="o">)</span>
<span class="k">done</span>
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="parallel">

<pre class="highlight"><code>stringtie_on_one_bam<span class="o">()</span> <span class="o">{</span>
    dx-download-all-inputs
    gunzip <span class="s2">"</span><span class="nv">$gtfgz_path</span><span class="s2">"</span>
    samtools sort -@ <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> -o sorted.bam
    mv sorted.bam <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span>
    mkdir -p ./out/assembled_transcripts_gtf
    stringtie <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> -o ./out/assembled_transcripts_gtf/<span class="s2">"</span><span class="nv">$bam_prefix</span><span class="s2">"</span>.assembled_transcripts.gtf <span class="se">\</span>
        -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> -G <span class="s2">"</span><span class="k">${</span><span class="nv">gtfgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span>

    dx-upload-all-outputs
<span class="o">}</span>
</code></pre>
  
  </div>
</div>



<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#single2" data-toggle="tab">single</a></li>
    <li><a href="#parallel2" data-toggle="tab">parallel</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="single2">

<pre class="highlight"><code>stringtie -G <span class="s2">"</span><span class="k">${</span><span class="nv">reference_gtfgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span> --merge <span class="s2">"</span><span class="k">${</span><span class="nv">assembled_transcripts_gtfs</span><span class="p">[@]</span><span class="k">}</span><span class="s2">"</span> <span class="se">\</span>
    -o ./out/merged_assembled_transcripts_gtf/merged.assembled_transcripts.gtf -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span>
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="parallel2">

<pre class="highlight"><code>stringtie_merge<span class="o">()</span> <span class="o">{</span>

    dx-download-all-inputs
    
    gunzip <span class="s2">"</span><span class="nv">$guide_gtfgz_path</span><span class="s2">"</span>
    mkdir -p ./out/merged_assembled_transcripts_gtf
    stringtie -G <span class="s2">"</span><span class="k">${</span><span class="nv">guide_gtfgz_path</span><span class="p">%.gz</span><span class="k">}</span><span class="s2">"</span> --merge <span class="s2">"</span><span class="k">${</span><span class="nv">assembled_transcripts_gtfs_path</span><span class="p">[@]</span><span class="k">}</span><span class="s2">"</span> <span class="se">\</span>
        -o ./out/merged_assembled_transcripts_gtf/merged.assembled_transcripts.gtf -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span>
    
    dx-upload-all-outputs
<span class="o">}</span>
</code></pre> 

 
  </div>
</div>



<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#single3" data-toggle="tab">single</a></li>
    <li><a href="#parallel3" data-toggle="tab">parallel</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="single3">

<pre class="highlight"><code><span class="k">for </span>idx <span class="k">in</span> <span class="s2">"</span><span class="k">${</span><span class="p">!bams_path[@]</span><span class="k">}</span><span class="s2">"</span>; <span class="k">do

    </span><span class="nv">bam_path</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">bams_path</span><span class="p">[</span><span class="nv">$idx</span><span class="p">]</span><span class="k">}</span><span class="s2">"</span>
    <span class="nv">bam_prefix</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">bams_prefix</span><span class="p">[</span><span class="nv">$idx</span><span class="p">]</span><span class="k">}</span><span class="s2">"</span>

    <span class="nv">ctab_output_path</span><span class="o">=</span>/home/dnanexus/stringtie_output/<span class="s2">"</span><span class="nv">$bam_prefix</span><span class="s2">"</span>.ctab
    ctab_output_paths+<span class="o">=(</span><span class="s2">"</span><span class="nv">$ctab_output_path</span><span class="s2">"</span><span class="o">)</span>
    stringtie -eB -G ./out/merged_assembled_transcripts_gtf/merged_assembled_transcripts.gtf <span class="se">\</span>
        <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> -o /home/dnanexus/stringtie_output/<span class="s2">"</span><span class="nv">$bam_prefix</span><span class="s2">"</span>.gtf
<span class="k">done</span> 
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="parallel3">

<pre class="highlight"><code>
stringtie_prepare_ballgown<span class="o">()</span> <span class="o">{</span>
    
    dx-download-all-inputs
    samtools sort -@ <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> -o sorted.bam
    mv sorted.bam <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span>

    stringtie -eB -G <span class="s2">"</span><span class="nv">$merged_assembled_transcripts_gtf_path</span><span class="s2">"</span> <span class="s2">"</span><span class="nv">$bam_path</span><span class="s2">"</span> <span class="se">\</span>
        -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> -o /home/dnanexus/stringtie_output/output.gtf

    mkdir -p ./out/ballgown_ctabs_targz
    <span class="nb">cd </span>stringtie_output <span class="o">&amp;&amp;</span> tar czvf ./out/ballgown_ctabs_targz/<span class="s2">"</span><span class="nv">$bam_prefix</span><span class="s2">"</span>.ballgown_ctabs.tar.gz <span class="k">*</span>.ctab <span class="o">&amp;&amp;</span> <span class="nb">cd</span> ..
    dx-upload-all-outputs
<span class="o">}</span>
</code></pre>

  
  </div>
</div>

Now that we have these individual functions that each are respsonsible for one piece of the StringTie
taks, we create a main function responsible for orchestrating the execution of the functions:

```shell
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
```

There are two DNAnexus features used in the `main` function. The first is the `dx-jobutil-new-job`
command. This is similar to `dx run`, but it runs a function on a separate worker as a "subjob" of
the current job. It also returns a job ID for that subjob that we can use to refer to its outputs.

The second feature is a "job-based object reference", which is the `$job_id:job_output_field` syntax
used in the `main` function. This allows us to refer to an output of a subjob before the subjob itself
has finished. By referring to these future values, we can set up all the jobs we need to run and their
dependencies immediately, so in this case the `main` function will complete very quickly, and all the
subjobs that it launches will evenutally supply its outputs.

## 9.2 Distributing in Python

The method for distributing tasks across workers in python is very similar to bash; the differences
are largely in syntax.

In python, subjobs are launched with `dxpy.new_dxjob`. Job-based object references are made with
`dxpy.DXJob.get_output_ref`. Finally, functions that will run on separate workers need to be
decorated with `dxpy.entry_point`.
