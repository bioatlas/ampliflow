#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'bioatlas/ampliflow': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'dada2filter.R': ['v_dada2filter.txt', r"dada2filter version (\S+)"]
    'dada2errmodels.R': ['v_dada2errmodels.txt', r"dada2errmodels version (\S+)"]
    'dada2cleanNmerge.R': ['v_dada2cleanNmerge.txt', r"dada2cleanNmerge version (\S+)"]
    'dada2bimeras.R': ['v_dada2bimeras.txt', r"dada2bimeras version (\S+)"]
    'dada2idseq.R': ['v_dada2idseq.txt', r"dada2idseq version (\S+)"]
    'dada2taxonomy.R': ['v_.txt', r"dada2taxonomy version (\S+)"]
}
results = OrderedDict()
results['bioatlas/ampliflow'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['dada2filter.R'] = '<span style="color:#999999;\">N/A</span>'
results['dada2errmodels.R'] = '<span style="color:#999999;\">N/A</span>'
results['dada2cleanNmerge.R'] = '<span style="color:#999999;\">N/A</span>'
results['dada2bimeras.R'] = '<span style="color:#999999;\">N/A</span>'
results['dada2idseq.R'] = '<span style="color:#999999;\">N/A</span>'
results['dada2taxonomy.R'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'bioatlas/ampliflow Software Versions'
section_href: 'https://github.com/bioatlas/ampliflow'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
