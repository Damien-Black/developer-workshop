---
title: Improving Execution
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_improving_execution.html
---

There are a number of ways to configure applets to run beyond the defaults used by the
platform.

## 6.0 Instance Type

If you do not specify an instance type for your applet, it will default to an instance
called "mem2_hdd2_x2". The user can override this, but it is better to specify a better
default value so the user does not have to worry about it.

This is done in the dxapp.json file:

```json
  ...
  "runSpec": {
    ...
    "systemRequirements": {"*": {"instanceType": "mem1_ssd1_x8"}},
  }
  ...
```

The list of DNAnexus instance types and the resources available for each one is available on
the DNAnexus [wiki](https://wiki.dnanexus.com/API-Specification-v1.0.0/Instance-Types).

### Selecting an instance type

Choosing the best default value for instance type is a bit of an art. There are a few
things to keep in mind:

- **Memory requirements**. You need to make sure that the worker has enough
  memory to accomodate the applet's highest memory usage.
- **Local storage requirements**. You have to make sure that there is enough
  disk space on the worker to fit input and output files for the applet.
- **Parallelization of multiple cores**. Generally, you would like to have your
  jobs finish faster, and you can often do that by increasing the number of
  CPUs that your applet uses. However, the improvement is not always linear.
  A program may only run, say, 20% faster if we double the number of CPUs.
- **Cost**. The cost of an instance increases with CPUs, memory, and storage.
  We want to avoid paying for resources that we are not going to use.

## 6.1 Selecting an OS Release

The default operating system use by DNAnexus workers is Ubuntu 12.04. However, Ubuntu 12.04
will reach end-of-life in April 2017, which is rather soon. This means that it will stop
receiving security updates and more generally will be increasingly out of date.

Ubuntu 12.04 will remain the default, but you can specify Ubuntu 14.04 in the dxapp.json file:

```json
  ...
  "runSpec": {
    ...
    "distribution": "Ubuntu",
    "release": "14.04",
  }
  ...
```

One caveat to using Ubuntu 14.04 is that not all instance typs on DNAnexus currently support it.
Specifically, any instance type with "hdd2" in its name will not work with Ubuntu 14.04.

But, unless you need one of those instance types, it is almost always better to switch to Ubuntu 14.04.

## 6.2 Setting a Timeout Policy

Sometimes an applet can enter a state where productive work has stopped, but there has been
no error and the script is not complete. For example, an executable can just hang rather than
fail. In that case, the applet will just keep running for 30 days or until you notice it and
terminate it. That can lead to an unexpectedly high cost for a failed job.

You can protect yourself against this by using a "timeoutPolicy" in the applet's dxapp.json:


```json
  ...
  "runSpec": {
    ...
    "timeoutPolicy": {"*": {"hours": 6, "minutes": 30}}
  }
  ...
```

Choose an timeout value that is longer than you reasonably expect the applet to run on a large input.

## 6.3 Versioning

Finally, you can indicate a version of your applet in the dxapp.json file. This has more of an impact
when you switch from applets to apps, but it is still a nice thing to keep track of in an applet:


```json
{
  "name": "hisat2",
   "version": "0.0.1",
  ...
```
