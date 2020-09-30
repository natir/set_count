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

### Get uniq kmer

```
jellyfish count -m 31 -s {number of kmer estimation} -t {threads} -o uniq_kmer.jf {reference}
jellyfish dump -L 0 -U 1 -c -o uniq_kmer.csv uniq_kmer.jf
```

### MPHF

#### Index

- k: size of kmer you want count (not upper than 32)
- index: path where index is write
- threads: number of threads set_count can use durring index building (not upper than 255)

```
set_count mphf_index {k} uniq_kmer.csv {index} {threads}
```

#### Count

- reads: a sequence file, in fasta or fastq compress in gzip or not, whose kmers you want to count
- counts: path where count is write

```
set_count mphf_count {index} {reads} {counts}
```

#### Dump

```
set_count mphf_dump {counts} {reference}
```

### MQF

#### Index

MQF_count count only kmer uniq in reference.

- k: size of kmer you want count (not upper than 32)
- number_of_kmer: number of line in uniq_kmer.csv multiply by 10
- out: prefix of output

```
set_count mqf_index {k} {number_of_kmer} uniq_kmer.csv {out}.mqf
```

#### Count

- reads: a sequence file, in fasta or fastq compress in gzip or not, whose kmers you want to count

```
cp {index}.mqf {count}.mqf
set_count mqf_count {k} {reads} {out}.mqf
```

#### Dump

```
set_count mqf_dump {k} {counts} {reference}
```
