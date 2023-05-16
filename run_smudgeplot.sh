#!/bin/bash
echo "/vol/biotools/bin/jellyfish dump -c -L 31 -U 10000 \
/vol/gf-yam/members/bpucker/20200214_genome_size_estimation/genome_size_data/fw_19/reads.jf \
| /vol/cluster-data/bpucker/bin/smudgeplot/exec/smudgeplot.py hetkmers -o /vol/gf-yam/members/bpucker/20200219_smudgeplot/kmer_pairs_19" \
| qsub -cwd -N smudge -l vf=200G -l arch=lx-amd64 -P rapresmabs -pe multislot 1 \
-o /vol/gf-yam/members/bpucker/20200219_smudgeplot/output.txt \
-e /vol/gf-yam/members/bpucker/20200219_smudgeplot/error.txt