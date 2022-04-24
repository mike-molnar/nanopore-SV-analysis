import sys
import csv
import argparse

parser = argparse.ArgumentParser( description='Calcualte the B-allele frequency for a longshot VCF file.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

in_fh = open(args.input, 'r')
out_fh = open(args.output, 'w')

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'locus', 'allele_freq'], delimiter='\t')

for record in csv_reader:
    frequency = record['allele_freq'].split(",")
    out_fh.write("%s\t%s\t%s\t%.2f\n" % (record['chromosome'], record['locus'], record['locus'], int(frequency[1])/(int(frequency[0]) + int(frequency[1]))))
