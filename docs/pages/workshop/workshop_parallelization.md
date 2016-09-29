---
title: Parallelization
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_parallelization.html
---

By "parallelization", we mean running processes in parallel within a single worker
to take advantage of multiple CPUs.

## 8.0 Setting Tool Parameters

The most obvious way to parallelize your applet is to make sure that tools that already
can take advantage of multiple CPUs are set to do so. For hisat2, there is a `-p <proc>`
parameter that tells it the number of processes to use. Generally, you want to have a tool
run using all the CPUs available on whatever instance it is running on, so rather than
hard coding the number of processes, you should use `nproc` or `multiprocessing.cpu_count()`
in bash and python, respectively.


<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#python" data-toggle="tab">python</a></li>
    <li><a href="#bash" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="python">
<pre class="highlight"><code><span class="n">hisat2_cmd_template</span> <span class="o">=</span> <span class="p">(</span><span class="s">"hisat2 --dta -p {nproc} -x {index_basename} -1 {mate1_fastq} "</span>
                        <span class="s">"-2 {mate2_fastq}"</span><span class="p">)</span>
<span class="n">hisat2_cmd</span> <span class="o">=</span> <span class="n">hisat2_cmd_template</span><span class="o">.</span><span class="n">format</span><span class="p">(</span>
    <span class="n">nproc</span><span class="o">=</span><span class="n">multiprocessing</span><span class="o">.</span><span class="n">cpu_count</span><span class="p">(),</span>
    <span class="n">index_basename</span><span class="o">=</span><span class="n">index_basename</span><span class="p">,</span>
    <span class="n">mate1_fastq</span><span class="o">=</span><span class="n">inputs_dict</span><span class="p">[</span><span class="s">"mate1_fastqgz_path"</span><span class="p">][</span><span class="mi">0</span><span class="p">][:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".gz"</span><span class="p">)],</span>
    <span class="n">mate2_fastq</span><span class="o">=</span><span class="n">inputs_dict</span><span class="p">[</span><span class="s">"mate2_fastqgz_path"</span><span class="p">][</span><span class="mi">0</span><span class="p">][:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".gz"</span><span class="p">)])</span>
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="bash">

<pre class="highlight"><code>hisat2 --dta -p <span class="s2">"</span><span class="k">$(</span>nproc<span class="k">)</span><span class="s2">"</span> -x <span class="nv">$index_basename</span> -1 <span class="s2">"</span><span class="nv">$mate1_fastq_path</span><span class="s2">"</span> -2 <span class="s2">"</span><span class="nv">$mate2_fastq_path</span><span class="s2">"</span> -S hisat2_output.sam
</code></pre>
  
  </div>
</div>

## 8.1 Parallelization in Bash

Commands can be run in parallel in bash by adding a `&` symbol at the end of the command. So if we had
this line in our script

```shell
( dx cat "$mate1_fastqgz" | gunzip > "$mate1_fastq_path" ) &
```

the download and decompression would begin, but the script would not wait for it to complete. Instead,
it will move on to the next line in the script. If we later need to be sure that this process is
finished before we continue, we can use the `wait` command along with the id of the process.

```shell
( dx cat "$mate1_fastqgz" | gunzip > "$mate1_fastq_path" ) &
download_pid=$!

...

wait "$download_pid"
```

Alternatively, we could use the `wait -n` command with will way for *any* background process to finish
and return its exit status.

{% include antipattern.html content="There are two forms of `wait` that are dangerous. The first is a `wait` with no process IDs. This will wait for all background jobs to finish, but it will ignore any errors. Similarly, `wait <pid1> <pid2> ...` with multiple process IDs will ignore errors in any process except the last one." %}

Putting this together, we could parallelize the download of the two FASTQ files at the beginning
of the hisat2 applet script:

```shell
    mate1_fastq_path="${mate1_fastqgz_path%.gz}"
    mkdir -p "$(dirname "$mate1_fastq_path")"
    ( dx cat "$mate1_fastqgz" | gunzip > "$mate1_fastq_path" ) &

    mate2_fastq_path="${mate2_fastqgz_path%.gz}"
    mkdir -p "$(dirname "$mate2_fastq_path")"
    ( dx cat "$mate2_fastqgz" | gunzip > "$mate2_fastq_path" ) &

    wait -n
    wait -n
```

This will download both files and extract them in parallel. Then when both are finished, the
script will continue.

{% include antipattern.html content="Do not try to run many downloads or uploads in parallel. `dx download` and `dx upload` each parallelize file transfer internally, so running many of them in parallel stresses the platform and can degrade performance." %}

Finally, though it is beyond the scope of this workshop, [GNU Parallel](https://www.gnu.org/software/parallel/) is a flexible tool
for running commands in parallel.

## 8.2 Parallelization in Python

The first way to parallelize commands in python is to use `subprocess.Popen`, whic will return immediately. Then, later you
use `subprocess.Popen.wait()` to wait for the command to complete. This is analogous to `&` and `wait` in bash. Remember
to check the return codes of the commands:

```python
gunzip_proc_1 = subprocess.Popen(["gunzip", inputs_dict["mate1_fastqgz_path"][0]])
gunzip_proc_2 = subprocess.Popen(["gunzip", inputs_dict["mate2_fastqgz_path"][0]])

gunzip_proc_1.wait()
gunzip_proc_2.wait()

if any([p.returncode != 0 for p in [gunzip_proc_1, gunzip_proc_2]]):
    raise dxpy.AppError("Gunzipping FASTQ files failed.")
```

The second way is to use the `multiprocessing` module. There are many ways to use this module, but the
most frequent in DNAnexus applets is to use a `multiprocessing.Pool`, which creates pool of processes
that can run functions in parallel.

```python
def unzip_file(file_path):
    try:
        subprocess.check_call(["gunzip", file_path])
    except Exception as exc:
        traceback.print_exc()
        raise Exception()

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastqgz, mate2_fastqgz):

    # First, download all the input files to local storage
    inputs_dict = dxpy.download_all_inputs()

    pool = multiprocessing.Pool() 
    pool.map(unzip_file, [inputs_dict["mate1_fastqgz_path"][0], inputs_dict["mate2_fastqgz_path"][0]])
```

In this code, the two fastq.gz paths are submitted to the `unzip_file1` function in parallel. 

{% include antipattern.html content="Do not run `subprocess.check_call` or `subprocess.check_output` using `multiprocessing` unless it is wrapped in a `try`-`except` block that will raise a bare `Exception`. Otherwise, if the command fails, the script will hang indefinitely." %}
