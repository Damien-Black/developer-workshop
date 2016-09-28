---
title: Organizing the Pipeline
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_organizing_the_pipeline.html
---

## 2.0 Pipeline Structure

Before we start any coding, we need to think about how to structure the pipeline
on the DNAnexus platform. The authors of the pipeline have provided a schematic
of its overall organization:

{% include image.html file="DE_pipeline.png" caption="" alt="Pipeline Schematic" %}

There are five different commands that are run within the pipeline, along with
some different handling of multiple files and samples. There are lots of ways
we could break this pipeline up into different applets on the platform. On one
extreme, we could write one applet that runs the whole pipeline. That would
mean that one DNAnexus worker would be responsible for running the whole pipeline.
On the other, we could have an applet for each command or even for command/input
combinations. That would distribute the pipeline into many potentially small
tasks run on different workers. When deciding how to structure a pipeline, keep
in mind the tradeoffs associated with creating more or fewer applets.

### Advantages to more applets
- **Better recoverability and restartability.** When a pipeline is broken into
  lots of different steps, it is easier to recover and restart and intermediate
  step that fails. On the other hand, if the pipeline is contained within a single
  applet, a failure can mean that the pipeline must be rerun from the beginning.

- **Better resource tuning.** When different steps of a pipeline are broken
  into different applets, each applet can be assigned a different instance type.
  For example, a step that needs lots of memory can get a high-memory instance,
  but a step that has lower memory requirements can be run on a cheaper instance
  with less memory.

- **Better modularity.** Generally, breaking a pipeline into different applets
  can make it easier for a developer to understand and implement each step, and
  it can make it easier to swap or update individual steps.

### Advantages to fewer applets

- **Lower file I/O.** Each applet runs on a separate worker, and workers each
  have independent local storage. So, having fewer applets can reduce file transfer
  between the platform and workers. For example if one step of a pipeline produces
  a file consumed by a second step, running both steps within the same applet will
  allow the second step to start immediately without having to download the output
  of the first step.

- **Lower overhead.** Each worker needs some time to get started. There can be
  some delay in requesting a worker from the infrastructure provider, and a worker
  has some initialization work to do before it can run the applet's code. So,
  having fewer applets reduces this execution overhead.

## 2.1 Dividing Our Pipeline into Applets

For our pipeline, we will implement three applets:

{% include image.html file="DE_pipeline_broken_up.png" caption="" alt="Pipeline Division" %}

So, one applet will be responsible for running HISAT2, one applet will run StringTie,
and one will run Ballgown. This lets us reduce the file I/O needed to pass the
alignments and assembled transcripts to the multiple StringTie commands. Also, since
Ballgown is an R application, it lets us keep R dependencies limited to one applet.

{% include links.html %}
