---
title: Welcome to the DNAnexus Developer Workshop
keywords: getting_started
sidebar: tutorial_sidebar
permalink: index.html
---
## 0.0 Welcome

Welcome to the DNAnexus Connect Developer Workshop!

This workshop will help you write tools on the DNAnexus platform more effectively. We are going to write applets that

- **Run efficiently**, reducing the cost and time needed to run analyses.
- **Are easy to debug**, letting developers understand and resolve issues
- **Use the scale of the cloud**, taking advantage of the DNAnexus platform's
  flexibility
- **Are easy to use**, reducing support and enabling collaboration 

## 0.1 Prerequisites

To get the most out of this workshop, you will need a couple things:

- A DNAnexus account so you can follow along and build and run applets.
- Some experience using the platform and writing applets. This is not strictly
  necessary, but we will elide some DNAnexus basics. However, there will be lots
  of code examples and a pre-built environment, so you can jump in wherever you
  feel comfortable.

## 0.2 Setting Up

1. Log in to the DNAnexus platform at platform.dnanexus.com, and tell me the
DNAnexus user name that you will be using, so I can share the workshop
project with you.
2. Create a project that will be the "sandbox" where you can build and test
applets.

The next steps depend on whether you are comfortable installing the SDK and
using it from the computer that you've brought with you to the workshop.

If you are comfortable doing that, then
1. Download the DNAnexus SDK if you do not already have it. You can get it
from [https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK](https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK).
Follow the instructions there to set up the SDK and verify that it is
working.
2. Clone the GitHub repo that contains the source files for the applets
that we will be working on: [https://github.com/dnanexus/developer-workshop](https://github.com/dnanexus/developer-workshop)

If you think you may have issues with installing and running the SDK, there
is an alternative. This may be the case if you know you have are using a
corporate Windows computer or you have an unusual distribution or Linux or
you are restricted from using new software.

1. Go to the DNAnexus Developer Workshop project through the DNAnexus
web UI.
2. Click on the "workshop-in-a-box" applet, and copy it to your new
sandbox project.
3. Run the applet with SSH access enabled.
4. Wait a minute or so after the applet starts, and check the log. There
you will see instructions to paste a URL into your browser. It will look
something like this:
https://123.456.789.00:22
5. You will then need to start your browser in a mode where it will allow
you to connect to port 22. This varies by browser:

*Google Chrome* - Run with `--explicity-allowed-ports=22`

*Firefox* - Follow the instructions [here](https://support.mozilla.org/en-US/questions/1083282)

6. When you go to the URL, your broswer will give you a warning about a certificate; choose to proceed
anyway. Then, you will see a command line that is running on a remote DNAnexus
worker.
7. Source the DNAnexus environment by running `source dx-toolkit environment`.
8. In the DNAnexus web UI, create an access token that has contribute access to the
sandbox project your created earlier.
9. Run `dx login --token <token>` to login.

Now you are ready to start the workshop.

{% include links.html %}
