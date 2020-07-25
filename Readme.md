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
