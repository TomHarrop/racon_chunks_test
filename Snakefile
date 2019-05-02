#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

racon_chunks = 'shub://TomHarrop/singularity-containers:racon-chunks'
py36 = 'shub://TomHarrop/singularity-containers:py3.6.8_biopython1.73'
py37 = 'shub://TomHarrop/singularity-containers:py3.7.3_biopython1.73'
py36_mod = 'shub://TomHarrop/singularity-containers:py3.6.8_biopython1.73_mod'
py37_mod = 'shub://TomHarrop/singularity-containers:py3.7.3_biopython1.73_mod'
py371 = 'shub://TomHarrop/singularity-containers:py3.7.1_biopython1.73'
py371_mod = 'shub://TomHarrop/singularity-containers:py3.7.1_biopython1.73_mod'


all_reads = 'all_reads.fq'

########

# MAIN #
########

# all_chunks = glob_wildcards('sam/chunk_{chunk}.sam').chunk
# some_chunks = [x for x in all_chunks if int(x) < 10]

fraction_to_map = 1
seed = 14


#########
# RULES #
#########

rule target:
    input:
        expand('py{py}{mod}/r1.idx',
               py=['36', '37', '371'],
               mod=['', '_mod'])
        # expand('fq/chunk_{chunk}.fq',
        #        chunk=some_chunks)

# rule repair_reads:
#     input:
#         r1 = 'fq/chunk_{chunk}_r1.fq',
#         r2 = 'fq/chunk_{chunk}_r2.fq'
#     output:
#         'fq/chunk_{chunk}.fq'
#     shell:
#         'repair.sh '
#         'in={input.r1} '
#         'in2={input.r2} '
#         'repair=t '
#         'out={output}'

# rule retrieve_reads:
#     input:
#         sam = 'sam/chunk_{chunk}.sam',
#         r1_idx = 'py37/r1.idx',
#         r2_idx = 'py37/r2.idx',
#         r1 = 'subset_r1.fq',
#         r2 = 'subset_r2.fq'
#     output:
#         r1 = temp('fq/chunk_{chunk}_r1.fq'),
#         r2 = temp('fq/chunk_{chunk}_r2.fq')
#     benchmark:
#         'fq/retrieve_{chunk}_benchmark.txt'
#     script:
#         'src/retrieve_reads.py'


rule index_reads_371_mod:
    input:
        'reads/r{r}.fq'
    output:
        'py371_mod/r{r}.idx'
    log:
        'py371_mod/r{r}.log'
    benchmark:
        'py371_mod/r{r}_benchmark.txt'
    singularity:
        py371_mod
    script:
        'src/index_reads.py'


rule index_reads_36_mod:
    input:
        'reads/r{r}.fq'
    output:
        'py36_mod/r{r}.idx'
    log:
        'py36_mod/r{r}.log'
    benchmark:
        'py36_mod/r{r}_benchmark.txt'
    singularity:
        py36_mod
    script:
        'src/index_reads.py'

rule index_reads_37_mod:
    input:
        'reads/r{r}.fq'
    output:
        'py37_mod/r{r}.idx'
    log:
        'py37_mod/r{r}.log'
    benchmark:
        'py37_mod/r{r}_benchmark.txt'
    singularity:
        py37_mod
    script:
        'src/index_reads.py'


rule index_reads_36:
    input:
        'reads/r{r}.fq'
    output:
        'py36/r{r}.idx'
    log:
        'py36/r{r}.log'
    benchmark:
        'py36/r{r}_benchmark.txt'
    singularity:
        py36
    script:
        'src/index_reads.py'

rule index_reads_37:
    input:
        'reads/r{r}.fq'
    output:
        'py37/r{r}.idx'
    log:
        'py37/r{r}.log'
    benchmark:
        'py37/r{r}_benchmark.txt'
    singularity:
        py37
    script:
        'src/index_reads.py'

rule index_reads_371:
    input:
        'reads/r{r}.fq'
    output:
        'py371/r{r}.idx'
    log:
        'py371/r{r}.log'
    benchmark:
        'py371/r{r}_benchmark.txt'
    singularity:
        py371
    script:
        'src/index_reads.py'


rule split_reads:
    input:
        all_reads
    output:
        r1 = 'reads/r1.fq',
        r2 = 'reads/r2.fq'
    params:
        fraction = fraction_to_map,
        seed = seed
    log:
        'reads/split.log'
    singularity:
        racon_chunks
    shell:
        'reformat.sh '
        'in={input} '
        'int=t '
        'verifyinterleaved=t '
        'samplerate={params.fraction} '
        'sampleseed={params.seed} '
        'out={output.r1} '
        'out2={output.r2} '
        '2> {log}'
