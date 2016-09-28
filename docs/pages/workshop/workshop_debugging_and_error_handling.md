---
title: Debugging and Error Handling
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_debugging_and_error_handling.html
---

## 5.0 Debug Hold 

The HISAT2 jobs that we just started should all complete successfully. However
if we use a different reference, `grch38_snp_tran.tar.gz`, we get an error:

{% include image.html file="failed_job_summary.png" caption="" alt="Failed Job Summary" %}

The error invites us to look at the log:

{% include image.html file="failed_job_log.png" caption="" alt="Failed Job Log" %}

None of that is terribly helpful, but we have not yet used any DNAnexus debugging or error
handling features. The first and most useful feature is "debug hold". This is a setting
used when running a job that tells the platform to keep the worker running for 3 days
if it encounters an error. This lets us SSH into the worker that was running the failed
job and figure out what happened.

We can run a job in debug hold by selecting it when running the job through the UI or
by adding this flag to the `dx run` call:

```shell
--debug-on All
```

Then, when the job enters the "debug hold" state, we can SSH into it:

```shell
$ dx ssh_config
```

and figure out what went wrong. In this case, we can look at the reference directory
that gets extracted and see that one of our assumptions was not correct:

{% include image.html file="debugon_error.png" caption="" alt="Debug Hold Error" %}

Then, we can fix our script to handle this case

### Bash

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old" data-toggle="tab">old</a></li>
    <li><a href="#fixed" data-toggle="tab">fixed</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old">

<pre class="highlight"><code>    <span class="c"># Get the index basename to use when calling hisat2</span>
    <span class="nv">index_filename</span><span class="o">=</span><span class="k">$(</span>dx describe --name <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span><span class="k">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index_filename</span><span class="p">%.tar.gz</span><span class="k">}</span>/genome
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="fixed">

<pre class="highlight"><code>    <span class="c"># Get the index basename to use when calling hisat2</span>
    <span class="nv">index_filename</span><span class="o">=</span><span class="k">$(</span>dx describe --name <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span><span class="k">)</span>
    <span class="nv">index_directory</span><span class="o">=</span><span class="k">${</span><span class="nv">index_filename</span><span class="p">%.tar.gz</span><span class="k">}</span>
    <span class="nv">index1_filename</span><span class="o">=(</span> <span class="nv">$index_directory</span>/<span class="k">*</span>.1.ht2 <span class="o">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index1_filename</span><span class="p">%.1.ht2</span><span class="k">}</span>
</code></pre>
  
  </div>
</div>

### Python

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#old2" data-toggle="tab">old</a></li>
    <li><a href="#fixed2" data-toggle="tab">fixed</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="old2">

<pre class="highlight"><code>    <span class="n">index_basename</span> <span class="o">=</span> <span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">dxpy</span><span class="o">.</span><span class="n">DXFile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">)</span><span class="o">.</span><span class="n">name</span><span class="p">[:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".tar.gz"</span><span class="p">)],</span>
                                  <span class="s">"genome"</span><span class="p">)</span>
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="fixed2">

<pre class="highlight"><code>    <span class="n">index_directory</span> <span class="o">=</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">DXFile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">)</span><span class="o">.</span><span class="n">name</span><span class="p">[:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".tar.gz"</span><span class="p">)]</span>
    <span class="n">index_basename</span> <span class="o">=</span> <span class="n">glob</span><span class="o">.</span><span class="n">glob</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">index_directory</span><span class="p">,</span> <span class="s">"*.1.ht2"</span><span class="p">))[</span><span class="mi">0</span><span class="p">][:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".1.ht2"</span><span class="p">)]</span>
</code></pre>
  
  </div>
</div>

## 5.1 Reporting Errors in Bash

The DNAnexus platform provides ways to gives users more helpful information when a job fails. 

### mark-section and mark-success
The first strategy is a pair of brief scripts called `mark-section` and `mark-success`.
These scripts record an error message that tells the user what the applet what doing
when it was failed. This replaces the "Please consult the job log" message. We can add
these to our bash script:

```shell
#!/bin/bash
set -x -o pipefail

main() {

    # Download the inputs to the worker's local storage
    mark-section "downloading input files"
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

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

    # Upload the resulting file and associate it with the aligned_sam output
    mark-section "uploading SAM results"
    uploaded_id=$(dx upload hisat2_output.sam --brief)
    dx-jobutil-add-output aligned_sam $uploaded_id
    mark-success
}
```

And of course, we need to add the `mark-section` and `mark-success` executables to the applet's
resources directory.

Now if the tar command fails, the user will see a message "Error while extracting reference tarball."

### `set -x -o pipefail`

We can also make bash logs more useful and make a bash script less likely to ignore an error. To do this, we add one line to the top of our bash script:

```shell
#!/bin/bash
set -x -o pipefail
```
The `-x` will have each line printed into the log as it is executed. This makes it easier to figure out exactly where in the applet script and error occurred.

The `-o pipefail` ensures that an error in the middle of a pipeline will not be ignored. In bash, you can pipe commands together, but when `pipefail` is not set, failures of commands in the pipeline will be ignored. For example if you had 

```shell
cmd1 | cmd2 | cmd3
```
and `cmd2` failed but `cmd3` succeeded, the script would continue without reporting an error. That is not usually what you want in a DNAnexus applet.

## 5.2 Reporting Errors in Python

Error reporting in python needs less special attention, but it can be customized by raising 
`dxpy.AppError` exceptions. This will change the error message that is presented to the user.
So, we could accomplish something similar to `mark-section` by wrapping parts of the script
in a `try`-`except` block:

```python
    try:
        proc = subprocess.Popen(["tar", "xf", "hisat2_index.tar.gz"])
        proc.wait()
    except Exception:
        raise dxpy.AppError("Error extracting reference tarball")
```

This gives a nicer message to the user. But, one issue to keep in mind is that in Python,
this can end up obscuring important debugging information in the log. So it is a good
idea to print the traceback if you catch a general exception:

```python
    try:
        proc = subprocess.Popen(["tar", "xf", "hisat2_index.tar.gz"])
        proc.wait()
    except Exception:
        traceback.print_exc()
        raise dxpy.AppError("Error extracting reference tarball")
```

### A common pitfall

In the code snippet above, there a few things that could go wrong with the tar command.
Maybe it runs out of disk space, or maybe there is a typo in the command. In those
cases, you would want the script to halt and report that the extraction failed. Unfortunately,
this code will not do that. It will just soldiering on, and you will likely see a misleading
error later in the script.

This happens because `subprocess.Popen.wait` does not raise an exception if the command
it runs has a nonzero exit code. If you want the script to fail if the command fails,
you either have to check the exit code directly, or you can use `subprocess.check_call`
or `subprocess.check_output` to run the command.

```python
    try:
        subprocess.check_call(["tar", "xf", "hisat2_index.tar.gz"])
    except Exception:
        traceback.print_exc()
        raise dxpy.AppError("Error extracting reference tarball")
```

## 5.3 Static Analysis Tools

Many bugs in applets can be discovered quickly without needing to wait to see it during an
actual job on the platform. Depending on the language you are using for your applet script,
you can use [shellcheck](https://github.com/koalaman/shellcheck) or
[pylint](https://www.pylint.org/) to find bugs in an automated way.

### Shellcheck

Shellcheck is an extremely useful tool when writing bash scripts on DNAnexus. There
are a lot of corner cases in bash and mistakes that are easy to make. Shellcheck finds
them and even has a nice [wiki](https://github.com/koalaman/shellcheck/wiki) that describes
the reasoning behind each issue it identifies.

Shellcheck can be apt-get intsalled in Ubuntu 14.04, and it available on the workshop-in-a-box.

```shell
$ shellcheck ~/applets/hisat2/bash/41_script.sh
```

### Pylint

Pylint is similar to shellcheck, except it is for python. It also tends to identify issues
that are more stylistic rather than something that would lead to a bug.

```shell
$ pylint ~/applets/hisat2/python/41_script.py
```
### dxapp_lint.py

Finally, you may want to ensure that the dxapp.json files you create follow a certain style
or meet certain requirements. There is a script at `resources/dxapp_lint.py` that checks a
dxapp.json file against a long list of rules. The script there is currently rather pedantic,
but you can examine the code and make adjustments relatively easily.

```shell
$ python ~/resources/dxapp_lint.py ~/applets/hisat2/bash/41_dxapp.json
```
