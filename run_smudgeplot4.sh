#!/bin/bash
echo "/vol/biotools/bin/jellyfish dump -c -L 28 -U 10000 \
/vol/gf-yam/members/bpucker/20200214_genome_size_estimation/genome_size_data/fw_25/reads.jf \
| /vol/cluster-data/bpucker/bin/smudgeplot/exec/smudgeplot.py hetkmers -o /vol/gf-yam/members/bpucker/20200219_smudgeplot/kmer_pairs_25" \
| qsub -cwd -N smudge -l vf=200G -l arch=lx-amd64 -P rapresmabs -pe multislot 1 \
-o /vol/gf-yam/members/bpucker/20200219_smudgeplot/output4.txt \
-e /vol/gf-yam/members/bpucker/20200219_smudgeplot/error4.txt
