"""Checks the styling and formatting of a dxapp.json file.

The goal of this script is to make it easier to maintain consistent
conventions when creating apps and applets. It runs different
checks on components of the dxapp.json file and reports potential issues.

As a general principle, this will err on the side of false positives. So,
some things it reports are going to be correct as they are.

Different checks are just implemented as functions with the "check"
decorator. The decorator accepts arguments that specify which elements
of the dxapp.json the function should be applied to. These are the valid
arguments:

    - dxapp: The entire parsed dxapp.json file as a dict
    - title: The title field of dxapp.json
    - summary: The summary field of dxapp.json
    - input: Objects in the inputSpec
    - output: Objects in the outputSpec
    - label: The label field of any input or output
    - help: The help field of any input or output

When a check function finds a problem, it should raise DXAppConventionError
with a helpful message that describes the problem and how to fix it.
"""

import argparse
import functools
import itertools
import json
import re
import string

REGISTERED_CHECKS = {}
def check(*check_types):
    """Little wrapper to keep track of what dxapp elements to apply checking
    functions to.
    """
    def wrap(func):
        for check_type in check_types:
            REGISTERED_CHECKS.setdefault(check_type, []).append(func)
        @functools.wraps(func)
        def wrapped_func(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapped_func
    return wrap

################################################################################
# Checks to run on dxapp.json
################################################################################

# Check that fields that make the applet easier to understand are present.
@check("dxapp")
def has_summary(dxapp_obj):
    """Applets should have a summary field that describes what it does."""
    if 'summary' not in dxapp_obj:
        raise DXAppConventionError(
            "The dxapp.json is missing a summary.")

@check("dxapp")
def has_title(dxapp_obj):
    """Applets should have a title so they look good in the UI."""
    if 'title' not in dxapp_obj:
        raise DXAppConventionError(
            "The dxapp.json is missing a title.")

@check("input", "output")
def has_label(io_obj):
    """Inputs and outputs should have a label."""
    if 'label' not in io_obj:
        raise DXAppConventionError(
            "The IO field '{n}' is missing a label."
            .format(n=io_obj['name']))

@check("input", "output")
def has_help(io_obj):
    """Inputs and outputs should have a help string."""
    if 'help' not in io_obj:
        raise DXAppConventionError(
            "The IO field '{n}' is missing a help message."
            .format(n=io_obj['name']))

# Check the formatting of customer-facing strings.
@check("label", "help", "summary", "title")
def first_letter_is_capitalized(string_):
    """The first character of user-facing strings should be capitalized.

    Or at least something should be capitalized, dbDNP, elPrep, for example."""
    if not any([s in string.ascii_uppercase for s in string_.split()[0]]):
        raise DXAppConventionError(
            "The first letter of '{s}' should be capitalized."
            .format(s=string_))

@check("help")
def ends_with_period(help_string):
    """Help strings should end with a period

    I guess. That seems to be the way things are done.
    """
    if help_string[-1] != '.':
        raise DXAppConventionError(
            "The string '{s}' should end in a period."
            .format(s=help_string))

@check("summary")
def summary_is_short(summary_string):
    """The summary should fit on one line when it shows up in the app store
    list, so it can't be too long.
    """
    if len(summary_string) > 120:
        raise DXAppConventionError(
            "The applet summary should be a one-liner, less than 120 characters not {c}."
            .format(c=len(summary_string)))

@check("input", "output")
def boolean_label_is_a_question(io_obj):
    """The labels of booleans should be yes or no questions."""
    if io_obj['class'] in ("boolean", "array:boolean"):
        if 'label' in io_obj and io_obj['label'][-1] != '?':
            raise DXAppConventionError(
                "The IO field {n} is a boolean, so its label "
                "should be a question with a yes or no answer."
                .format(n=io_obj['name']))

@check("input", "output")
def optional_in_help(io_obj):
    """If an IO field is optional, its help message should have "(Optional)" at
    the start.
    """
    if io_obj.get("optional"):
        if 'help' in io_obj:
            if not io_obj['help'].startswith("(Optional)"):
                raise DXAppConventionError(
                    "The optional IO field '{n}' should have a help message that "
                    "that starts with '(Optional)'."
                    .format(n=io_obj['name']))

@check("label", "summary", "help", "title")
def correct_spellings(check_string):
    """Use the canonical spelling/capitalization of common bioinformatics
    terms.
    """
    correct_spellings_ = {
        "bwa": "BWA",
        "bwa.mem": "BWA-MEM",
        "bwa.sw": "BWA-SW",
        "bwa.backtrack": "BWA-backtrack",
        "gatk": "GATK",
        "gzip": "gzip",
        "lofreq": "LoFreq",
        "fastq": "FASTQ",
        "fasta": "FASTA",
        "bam": "BAM",
        "sam": "SAM",
        "cram": "CRAM",
        "gtf": "GTF",
        "illumina": "Illumina",
        "pacbio": "PacBio",
        "dnanexus": "DNAnexus",
        "freebayes": "FreeBayes",
        "sambamba": "Sambamba",
        "snp": "SNP",
        "dbsnp": "dbSNP",
        "snv": "SNV",
        "grch37": "GRCh37",
        "grch38": "GRCh38",
        "hg19": "hg19",
        "vcf": "VCF"}

    for term, correct_spelling in correct_spellings_.iteritems():
        term_as_single_word = r'(?<=[\s^])' + term + r'(?=[\s\.\?$])'
        for matching_instance in re.findall(term_as_single_word, check_string):
            if matching_instance != correct_spelling:
                raise DXAppConventionError(
                    "The string '{s}' in '{u}' should be written as '{t}'."
                    .format(s=matching_instance, t=correct_spelling,
                            u=check_string))

# Check that some reliability-improving dxapp elements are included.
@check("dxapp")
def specifies_instance_type(dxapp_obj):
    """Don't rely on the default instance type. It's rarely the best choice,
    and when it is, it's better to be explicit about it.
    """
    run_spec = dxapp_obj['runSpec']
    try:
        has_instance_type = 'instanceType' in run_spec['systemRequirements']['main']
    except KeyError:
        has_instance_type = False
    if not has_instance_type:
        try:
            has_instance_type = 'instanceType' in run_spec['systemRequirements']['*']
        except KeyError:
            has_instance_type = False
    if not has_instance_type:
        raise DXAppConventionError(
            "The dxapp.json should specify an instance type rather than using the default.")

@check("dxapp")
def specifies_timeout(dxapp_obj):
    """Applets should have a timeout policy so they halt when it's obvious that something
    has gone wrong rather than run up a large AWS bill.
    """
    if 'timeoutPolicy' not in dxapp_obj['runSpec']:
        raise DXAppConventionError(
            "The dxapp.json should specify a timeoutPolicy.")

@check("dxapp")
def only_apt_execdepends(dxapp_obj):
    """For applets that need to usually work, non-apt dependencies are a bad idea."""
    if 'execDepends' in dxapp_obj['runSpec']:
        for dependency in dxapp_obj['runSpec']['execDepends']:
            if 'package_manager' in dependency and dependency['package_manager'] != 'apt':
                raise DXAppConventionError(
                    "Don't use non-apt dependencies, connections to other package "
                    "managers sporadically fail. Replace with an apt package or "
                    "bundledDepends.")

################################################################################
# End of checks
################################################################################

class DXAppConventionError(Exception):
    """Exception that checks should raise when a problem is found."""
    pass

def run_check(check_func, obj):
    """Apply a check function, catching and printing the convention error."""
    try:
        check_func(obj)
    except DXAppConventionError as exc:
        print exc.message

def run_checks(dxapp_obj):
    """Apply the registered checks to the appropriate elements of
    the parsed dxapp.json file.
    """
    for check_ in REGISTERED_CHECKS['dxapp']:
        run_check(check_, dxapp_obj)

    for check_ in REGISTERED_CHECKS['title']:
        if 'title' in dxapp_obj:
            run_check(check_, dxapp_obj['title'])

    for check_ in REGISTERED_CHECKS['summary']:
        if 'summary' in dxapp_obj:
            run_check(check_, dxapp_obj['summary'])

    for check_ in REGISTERED_CHECKS['input']:
        for input_ in dxapp_obj.get('inputSpec', []):
            run_check(check_, input_)

    for check_ in REGISTERED_CHECKS['output']:
        for output in dxapp_obj.get('outputSpec', []):
            run_check(check_, output)

    for check_ in REGISTERED_CHECKS['label']:
        for io_obj in itertools.chain(dxapp_obj.get("inputSpec", []),
                                      dxapp_obj.get("outputSpec", [])):
            if 'label' in io_obj:
                run_check(check_, io_obj['label'])

    for check_ in REGISTERED_CHECKS['help']:
        for io_obj in itertools.chain(dxapp_obj.get("inputSpec", []),
                                      dxapp_obj.get("outputSpec", [])):
            if 'help' in io_obj:
                run_check(check_, io_obj['help'])

def get_parser():
    """Return an argument parser for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("dxapp_json")
    return parser

def main():
    """Entry point."""
    parser = get_parser()
    args = parser.parse_args()

    with open(args.dxapp_json) as dxapp:
        dxapp_obj = json.load(dxapp)

    run_checks(dxapp_obj)

if __name__ == '__main__':
    main()
