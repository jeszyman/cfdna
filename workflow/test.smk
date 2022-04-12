import pandas as pd
import numpy as np

samples = pd.read_table("/home/jeszyman/repos/mpnst-preprocessing/test/inputs/samples.tsv")

test = dict(zip(samples['old_name'], samples['new_name']))

test['mpnst1']

test.keys()
test.values()


rule all:
    input:
        expand('test/symlink/{value}_R1.fastq.gz', value=test.values())

rule my_job:
    output:
        "test/symlink/{f}_R1.fastq.gz"
    params:
        key=lambda wcs: test[wcs.f]
    shell:
        """
        ln -s {params.key} {output}
        """
