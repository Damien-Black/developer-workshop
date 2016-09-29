---
title: Handling Complex Dependencies
keywords: getting_started
sidebar: tutorial_sidebar
permalink: workshop_handling_complex_dependencies.html
---

## 10.0 Ballgown's Dependencies

The last applet we will work on is for Ballgown. Ballgown is an R package that provides a
lot of flexibility in creating figures and other differential expression analyses. The
R scripts will need to begin with

```
library(ballgown)
```

So, we will need to install Ballgown before we can run any of the R scripts we may want to
run. The Ballgown instructions tell us to do this via bioconductor, which can be done
very simply:

```
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
```

One solution for creating our Ballgown applet would be to have an R script that installs
the package from bioconductor at the beginning of the applet script. While this would
often work, it would introduce a dependency on an external service that may frequently
be unavailable. This is not limited to bioconductor packages; packages installed via
pip, cran, cpan, etc. will all have some risk of failing due to a transient inability
to download packages.

Additionally, dependencies can take a long time to install. They may require lengthy
compilation or data transfer.

## 10.1 Ballgown Asset

We can solve this issues by using `assets` on DNAnexus. There are two steps to
using an asset:

1. Creating the asset using `dx build_asset`.
2. Referring to the build asset using `assetDepends` in the `dxapp.json` of the 
applet.

The asset directory contains a file called dxasset.json that is similar in spirit
to a dxapp.json. It also optionally has a Makefile that contains instructions for
dependencies that cannot be installed via a package manager. For the Ballgown asset,
the dxasset.json is this:

```json
{
  "name": "ballgown_asset",
  "title": "Ballgown",
  "description": "Asset containing the dependencies for Ballgown.",
  "version": "0.0.2",
  "distribution": "Ubuntu",
  "release": "14.04"
}
```

And the Makefile contains code for upgrading R and then installing Ballgown
with bioconductor:

```makefile
all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo apt-get --yes update
	sudo apt-get remove --auto-remove --yes r-base
	sudo apt-get install --yes --force-yes r-base
	sudo apt-get install --yes --force-yes r-base-dev
	R -e 'source("http://bioconductor.org/biocLite.R"); biocLite(); biocLite("ballgown")'
```

When we run `dx build_asset`, this will create a record which we can use the dxapp.json
file for the actual Ballgown applet:

```json
  "runSpec": {
    "file": "src/script.sh",
    "interpreter": "bash",
    "systemRequirements": {"*": {"instanceType": "mem1_ssd1_x4"}},
    "distribution": "Ubuntu",
    "release": "14.04",
    "timeoutPolicy": {"*": {"hours": 6, "minutes": 30}},
    "assetDepends": [{"id": "record-BzbY2z80YZX3Xx1Q1kKZ808V"}]
  }
```
