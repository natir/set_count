#!/usr/bin/bash 

set -x

jellyfish count -m 31 -s 5000000 -o uniq_kmer.jf reference.fasta
jellyfish dump -L 0 -U 1 -c -o uniq_kmer.csv uniq_kmer.jf

\time -v ./bin/set_count mphf_index 31 uniq_kmer.csv mphf.index 1
\time -v ./bin/set_count mphf_count mphf.index reads.fasta mphf.counts
\time -v ./bin/set_count mphf_dump mphf.counts reference.fasta > mphf.dump

\time -v ./bin/set_count mqf_index 31 $(($(wc -l uniq_kmer.csv | cut -d$' ' -f1) * 2)) uniq_kmer.csv mqf.index
cp mqf.index mqf.counts
\time -v ./bin/set_count mqf_count 31 reads.fasta mqf.counts
\time -v ./bin/set_count mqf_dump 31 mqf.counts reference.fasta > mqf.dump

\time -v jellyfish count -m 31 -s 20000000 -o jellyfish.counts reads.fasta
