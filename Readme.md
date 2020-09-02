# set_count

Count kmer present of a specific set of kmer

## Building

```
git clone --recursive https://github.com/natir/set_count.git
cd set_count
make
```

`set_count` binary is present in bin directory `set_count.a` library is present in lib directory

## Usage

### Index

- k: size of kmer you want count (not upper than 32)
- kmers: a sequence file, in fasta or fastq compress in gzip or not, contains kmer you want count, sequence short than k is ignore
- index: path where index is write
- threads: number of threads set_count can use durring index building (not upper than 255)

```
set_count index {k} {kmers} {index} {threads}
```

### Count

- reads: a sequence file, in fasta or fastq compress in gzip or not, whose kmers you want to count
- counts: path where count is write

```
set_count count {index} {reads} {counts}
```

### Dump

```
set_count dump {counts}
```

Count is write in stdout in csv format

### MQF_index

MQF_count count only kmer uniq in reference.

- k: size of kmer you want count (not upper than 32)
- threads: number of threads set_count can use durring index building (not upper than 255)
- out: prefix of output
- reads: a sequence file, in fasta or fastq compress in gzip or not, whose kmers you want to count

```
jellyfish count -m {k} -s {number of kmer estimation} -t {threads} -o {out}.jf reference.fasta
jellyfish dump -L 0 -U 1 -c -o {out}.csv {out}.jf
set_count mqf_count {k} {threads} {out}.csv {reads} {out}.mqf
```

### MQF_count

```
cp {index}.mqf {count}.mqf
set_count mqf_count {k} {reads} {count}.mqf
```
