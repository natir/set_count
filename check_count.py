#!/usr/bin/env python

import csv

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))

def main(first, second):
    first_count = dict()

    with open(first) as fh:
        reader = csv.reader(fh)
        for row in reader:
            first_count[row[0]] = int(row[1])

    with open(second) as fh:
        reader = csv.reader(fh)
        for row in reader:
            key = ""
            if row[0] not in first_count:
                if rev_comp(row[0]) not in first_count:
                    print(f"kmer {row[0]} {rev_comp(row[0])} not in first")
                    continue
                else:
                    key = rev_comp(row[0])
            else:
                key = row[0]

            if first_count[key] != int(row[1]):
                print(f"kmer {row[0]} is different {first_count[key]} {row[1]}")

if __name__ == '__main__':
    import sys

    main(sys.argv[1], sys.argv[2])
