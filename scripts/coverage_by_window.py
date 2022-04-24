import sys
import csv
import argparse

parser = argparse.ArgumentParser( description='Calculate mean coverage levels in a sorted BED file by window size.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-w', '--window', type=int, required=False, default=500000)
parser.add_argument('-c', '--coverage', type=float, required=True)
parser.add_argument('-p', '--ploidy', type=float, required=False, default=2.0)
args = parser.parse_args()

in_fh = open(args.input, 'r')
out_fh = open(args.output, 'w')

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'locus', 'depth'], delimiter='\t')

region_start = 0
region_end = args.window
total_coverage = 0
total_locus = 0

for record in csv_reader:
    if int(record['locus']) > region_end:
        while int(record['locus']) > region_end:
            if total_locus == 0:
                out_fh.write("%s\t%d\t%d\t%d\t%.2f\t%.2f\n" % (record['chromosome'], region_start, region_end, 0, 0.0,0.0))
            else:
                mean_coverage = total_coverage/total_locus
                copy_number = float((mean_coverage / args.coverage) * args.ploidy)
                normalized = float(mean_coverage/(args.coverage * args.ploidy))
                if normalized > 1.0:
                    normalized = 1.0
                out_fh.write("%s\t%d\t%d\t%d\t%.2f\t%.2f\n" % (record['chromosome'], region_start, region_end, mean_coverage, normalized, copy_number))
            total_coverage=0
            total_locus = 0
            region_start = region_end
            region_end = region_end + args.window
        total_coverage = int(record['depth'])
        total_locus = 1
    else:
        total_coverage = total_coverage + int(record['depth'])
        total_locus = total_locus + 1

