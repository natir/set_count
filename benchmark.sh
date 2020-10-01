#!/usr/bin/bash 

set -x

mkdir -p benchmark

jellyfish count -m 31 -s 5000000 -o benchmark/uniq_kmer.jf dataset/reference.fasta
jellyfish dump -L 0 -U 1 -c -o benchmark/uniq_kmer.csv benchmark/uniq_kmer.jf

\time -v ./bin/set_count mphf_index 31 benchmark/uniq_kmer.csv benchmark/mphf.index 1
\time -v ./bin/set_count mphf_count benchmark/mphf.index dataset/reads.fasta benchmark/mphf.counts
\time -v ./bin/set_count mphf_dump benchmark/mphf.counts dataset/reference.fasta > benchmark/mphf.dump

\time -v ./bin/set_count mqf_index 31 $(($(wc -l benchmark/uniq_kmer.csv | cut -d$' ' -f1) * 2)) benchmark/uniq_kmer.csv benchmark/mqf.index
cp benchmark/mqf.index benchmark/mqf.counts
\time -v ./bin/set_count mqf_count 31 dataset/reads.fasta benchmark/mqf.counts
\time -v ./bin/set_count mqf_dump 31 benchmark/mqf.counts dataset/reference.fasta > benchmark/mqf.dump

\time -v jellyfish count -m 31 -s 20000000 -o benchmark/jellyfish.counts dataset/reads.fasta
