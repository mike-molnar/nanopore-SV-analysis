import sys
import csv
import argparse
import pybedtools
from pybedtools import BedTool
from pybedtools import genome_registry
from pybedtools.featurefuncs import extend_fields

# define a function to process the copy number vaiants
def process_cnv(feature):
    # get the copy number information
    coverage = float(feature[3])
    copy_number = float((coverage / args.mean_coverage) * args.ploidy)
    length = int(feature[2]) - int(feature[1])
    
    # check for gains and/or losses
    gain = 0
    loss = 0
    cnv_list = feature[4].split(",")
    for current_cnv in cnv_list:
        if float((float(current_cnv) / args.mean_coverage) * args.ploidy) >= args.high_copy_threshold:
            gain += 1
        else:
            loss += 1
            
    if gain > 0 and loss == 0:
        gain_loss = "gain"
    if gain == 0 and loss > 0:
        gain_loss = "loss"
    if gain > 0 and loss > 0:
        gain_loss = "gains and losses"
    
    # the feature has to be expanded before adding the data
    feature = extend_fields(feature, 6)
    feature[3] = str(length)
    feature[4] = "{:.2f}".format(copy_number)
    feature[5] = str(gain_loss)
    
    return feature

# define a function to swap the left and right locations in a BEDPE file
def swap_left_and_right(feature):
    # save the rigth side in variables
    right_chr = feature[3]
    right_start = feature[4]
    right_end = feature[5]
    
    # move the left side to the right
    feature[3] = feature[0]
    feature[4] = feature[1]
    feature[5] = feature[2]
    
    # copy the right side to the left
    feature[0] = right_chr
    feature[1] = int(right_start)
    feature[2] = int(right_end)
    
    return feature
    
# parse the arguements from the user
parser = argparse.ArgumentParser(description='Filter structural variant calls from BED and BEDPE files.')
# define input BED files
parser.add_argument('-ins', '--insertions', type=str, required=False)
parser.add_argument('-del', '--deletions', type=str, required=False)
parser.add_argument('-inv', '--invertions', type=str, required=False)
parser.add_argument('-dup', '--duplications', type=str, required=False)
parser.add_argument('-trans', '--translocations', type=str, required=False)
parser.add_argument('-depth', '--coverage-depth', type=str, required=False)
parser.add_argument('-split', '--split-alignments', type=str, required=False)
parser.add_argument('-low_map', '--low-mapping-regions', type=str, required=False)
parser.add_argument('-gaps', '--genome-gaps', type=str, required=True)
# define output BED files
parser.add_argument('-ins_out', '--insertions-output', type=str, required=False)
parser.add_argument('-del_out', '--deletions-output', type=str, required=False)
parser.add_argument('-inv_out', '--invertions-output', type=str, required=False)
parser.add_argument('-dup_out', '--duplications-output', type=str, required=False)
parser.add_argument('-trans_out', '--translocations-output', type=str, required=False)
parser.add_argument('-cnv_out', '--copy-number-variant-output', type=str, required=False)
# define variables for filtering
parser.add_argument('-cov', '--mean-coverage', type=float, required=True)
parser.add_argument('-p', '--ploidy', type=int, required=False, default=2)
parser.add_argument('-hc', '--high-copy-threshold', type=float, required=False, default=2.7)
parser.add_argument('-lc', '--low-copy-threshold', type=float, required=False, default=1.3)
parser.add_argument('-min_len', '--min-length', type=float, required=False, default=50)
parser.add_argument('-min_cov_indel', '--min-coverage-indels', type=float, required=False)
parser.add_argument('-min_cov_trans', '--min-coverage-translocations', type=float, required=False)
parser.add_argument('-min_call', '--min-calls', type=int, required=False, default=3)
parser.add_argument('-slop', '--slop-pct', type=float, required=False, default=0.5)
args = parser.parse_args()

# set the coverage thresholds if they were not provided
if args.min_coverage_indels:
    min_cov_indels = args.min_coverage_indels
else:
    min_cov_indels = args.mean_coverage * 0.10
    
if args.min_coverage_translocations:
    min_cov_translocations = args.min_coverage_translocations
else:
    min_cov_translocations = args.mean_coverage * 0.25

# store BED files into a BedTool object
if args.coverage_depth:
    in_depth = BedTool(args.coverage_depth)
if args.insertions:
    in_ins = BedTool(args.insertions)
if args.deletions:
    in_del = BedTool(args.deletions)
if args.invertions:
    in_inv = BedTool(args.invertions)
if args.duplications:
    in_dup = BedTool(args.duplications)
if args.translocations:
    in_trans = BedTool(args.translocations)
if args.genome_gaps:
    in_gaps = BedTool(args.genome_gaps)
if args.split_alignments:
    in_split = BedTool(args.split_alignments)
if args.low_mapping_regions:
    in_low_map = BedTool(args.low_mapping_regions)

# create a BED file of filtered insertions
if args.insertions and args.insertions_output:
    in_ins.filter(lambda x: (float(x[4]) >= args.min_length and float(x[5]) >= min_cov_indels and float(x[6]) >= args.min_calls))\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .saveas(args.insertions_output)

# create a BED file of filtered deletions
if args.deletions and args.deletions_output:
    in_del.filter(lambda x: (float(x[4]) >= args.min_length and float(x[5]) >= min_cov_indels and float(x[6]) >= args.min_calls))\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .saveas(args.deletions_output)

# create a BED file of filtered invertions
if args.invertions and args.invertions_output:
    in_inv.filter(lambda x: (float(x[4]) >= args.min_length and float(x[5]) >= min_cov_indels and float(x[6]) >= args.min_calls))\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .saveas(args.invertions_output)

# create a BED file of filtered duplications
if args.duplications and args.duplications_output:
    in_dup.filter(lambda x: (float(x[4]) >= args.min_length and float(x[5]) >= min_cov_indels and float(x[6]) >= args.min_calls))\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .saveas(args.duplications_output)

# create a BED file of filtered copy number variants
if args.coverage_depth and args.copy_number_variant_output:
    filtered_cnv = in_depth\
    .filter(lambda x: ((float(x[3])/args.mean_coverage)*args.ploidy) <= args.low_copy_threshold or ((float(x[3])/args.mean_coverage)*args.ploidy) >= args.high_copy_threshold)\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .merge(d=1000000, c=(4,4), o=("mean","collapse"))
    
    # process the estimated copy number variants
    cnv_final = filtered_cnv.each(process_cnv)
    cnv_final.saveas(args.copy_number_variant_output)

# create a BED file of filtered translocations
if args.translocations and args.translocations_output and args.insertions and args.split_alignments:
    # create a BedTool of insertions larger than 200bp for filtering
    large_insertions = in_ins\
    .filter(lambda x: (float(x[4]) >= 200 and float(x[5]) >= min_cov_indels and float(x[6]) >= args.min_calls))

    # filter the left side of the BEDPE
    filtered_left_side = in_trans\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_split.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(large_insertions.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)
    
    # switch the right side to the left in order to filter the other side of the BEDPE
    right_side = filtered_left_side.each(swap_left_and_right)
    
    filtered_right_side = right_side\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_low_map.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(in_split.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\
    .intersect(large_insertions.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)
    
    # swap the left and right back to the original format
    trans_filtered = filtered_right_side.each(swap_left_and_right)
    trans_merged = trans_filtered\
    .filter(lambda x: (float(x[7]) >= min_cov_translocations))\
    .merge(d=100, c=(4,5,6,7,8), o=("distinct","min","max","count","mean"))\
    .filter(lambda x: (int(x[6]) >= args.min_calls))\
    .saveas(args.translocations_output)
