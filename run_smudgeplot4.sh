#!/bin/bash
echo "/vol/biotools/bin/jellyfish dump -c -L 28 -U 10000 \
/./20200214_genome_size_estimation/genome_size_data/fw_25/reads.jf \
| /./bin/smudgeplot/exec/smudgeplot.py hetkmers -o /./20200219_smudgeplot/kmer_pairs_25" \
| qsub -cwd -N smudge -l vf=200G -l arch=lx-amd64 -P xxx -pe multislot 1 \
-o /./20200219_smudgeplot/output4.txt \
-e /./20200219_smudgeplot/error4.txt
