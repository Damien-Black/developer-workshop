---
title: Getting Started
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_getting_started.html
---

## 1.0 Workshop Structure

We are going to illustrate applet writing concepts and best practices by
writing a set of applets that runs a real pipeline: the HISAT2, StringTie,
Ballgown RNA-seq quantification pipeline published in
[Nature Protocols](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html) last month.

{% include image.html file="nature_protocols_abstract.png" caption="" alt="Nature Protocols Abstract" %}

The pipeline comprises three tools:
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) for read alignment,
[StringTie](http://www.ccb.jhu.edu/software/stringtie/index.shtml) for
transcript assembly and quantification, and
[Ballgown](https://github.com/alyssafrazee/ballgown) for differential
expression analysis. The authors have provided a schematic describing how the
different tools fit together within their recommended pipeline:

While this particular pipeline may not be of interest to you, it is relatively
complex, which will give us the opportunity to cover quite a few concepts that
you can apply to tools that you do want to implement. It also lets us show how
you can start from a complicated description of a pipeline and end with a tool
that anyone can run in a few clicks.

## 1.1 Applet Source Code

We are going to start from a basic and suboptimal HISAT2 applet and keep
adding features. The code associated with each iteration is available
in the `applets` directory of the GitHub repo or in the `~/applets` directory
of the workshop-in-a-box applet. Each file is prefixed with a couple numbers
that indicate with part of the workshop the file corresponds to. For example,
if there were a bash script for this section, it would be called `11_script.sh`
for section 1.1. You can use these intermediate files to save yourself some typing,
to catch up or skip ahead, and to compare code between different iterations of
the applet.

## 1.2 N.B - Apps vs. Applets

A quick note for more experienced DNAnexus users. DNAnexus has several types of
executables: applets, apps, and workflows. We are going to build a workflow of
applets, and we will not discuss apps very much. Apps are, in some sense, a
wrapper around an underlying applet. And we are going to focus on the code
within the applet.

{% include links.html %}
