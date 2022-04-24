import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def is_number(n):
    try:
        float(n)
    except ValueError:
        return False
    return True

def to_sam(query_name, strand, target_name, position, qual, cigar):
    if(strand == '+'):
        bits = 0
    else:
        bits = 16
    print(query_name+"\t"+str(bits)+"\t"+target_name+"\t"+str(position)+"\t"+str(qual)+"\t"+cigar+"\t*\t0\t0\t*\t*")

# Removes soft and hard clips off the end of a cigar string
def clean_cigar(cigar):
    seperator = ''
    cgs = re.findall('[0-9]*[A-Z]', cigar)
    start = False
    end = False
    if cgs[0].endswith('H') or cgs[0].endswith('S'):
        start = True
    if cgs[len(cgs)-1].endswith('H') or cgs[len(cgs)-1].endswith('S'):
        end = True
    if start and end:
        return seperator.join(cgs[1:len(cgs)-1])
    elif start:
        return seperator.join(cgs[1:])
    elif end:
        return seperator.join(cgs[:len(cgs)-1])
    else:
        return seperator.join(cgs)

# Checks for an insertion larger than the given distance in a bam record
def check_insertion_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    if record.mapping_quality < args.minimum_mapping_qual:
        return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('I'):
            if int(cg[:cg.find("I")]) >= int(distance):
                return True
                break;
    return False

# Checks for an insertion or hard clip larger than the given distance in a bam record
# Needed as hard clips don't appear in the seq in a record. Need to get them from the fastq
def check_insertion_hard_clip_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    #if record.mapping_quality < args.minimum_mapping_qual:
    #    return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('H'):
            if int(cg[:cg.find("H")]) > 0:
                return True
                break;
    return False

# Checks for an insertion, hardclip or softclip larger than the given distance in a bam record
def check_insertion_or_soft_clip_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    #if record.mapping_quality < args.minimum_mapping_qual:
    #    return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('I'):
            if int(cg[:cg.find("I")]) >= int(distance):
                return True
        if cg.endswith('S'):
            if int(cg[:cg.find("S")]) >= int(distance):
                return True
        if cg.endswith('H'):
            if int(cg[:cg.find("H")]) >= int(distance):
                return True
    return False

# Gets the total read length based on the cigar string
def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
    return count

def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0

def get_start_end_string(cigar1, cigar2, read_id, ref1_start, ref1_end, ref2_start, ref2_end):
    # Assume the cigars are oriented such the alignment is cigar1 - indel - cigar2
    # Remove any hard or soft clipped bases from the start of the read. This is our starting position
    # Next iterate over a cleaned cigar string to get how many bases in first half on read, add the indel size and then iterate over second cleaned cigar string
    # now we know length of read to traverse. Gen ends pos by adding length to start pos
    cg1_start = 0
    cgs = re.findall('[0-9]*[A-Z]', cigar1)
    if cgs[0].endswith('S'):
        cg1_start += int(cgs[0][:cgs[0].find("S")])
    elif cgs[0].endswith('H'):
        cg1_start += int(cgs[0][:cgs[0].find("H")])
    cg1_cleaned = clean_cigar(cigar1)
    cg1_len = get_read_length(cg1_cleaned)
    cg1_end = cg1_start + cg1_len
    cg2_start = 0
    cgs = re.findall('[0-9]*[A-Z]', cigar2)
    if cgs[0].endswith('S'):
        cg2_start += int(cgs[0][:cgs[0].find("S")])
    elif cgs[0].endswith('H'):
        cg2_start += int(cgs[0][:cgs[0].find("H")])
    cg2_cleaned = clean_cigar(cigar2)
    cg2_len = get_read_length(cg2_cleaned)
    cg2_end = cg2_start + cg2_len
    if (cg1_start >= cg2_end) or ((cg2_start - cg1_end) < 0):
        return "0_0_0"
    if ref2_start - ref1_end == 0:
        # no issues
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(cg2_start - cg1_end) + "I" + cg2_cleaned
    elif ref2_start - ref1_end < 0:
        n_count1 = abs(ref2_start - ref1_end)/2
        n_count2 = n_count1
        if abs(ref2_start - ref1_end) % 2 == 1:
            n_count2 += 1
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(n_count2) + "N" + str((cg2_start - cg1_end)) + "I" + str(n_count1) + "N" + cg2_cleaned
    else:
        # Need to mofidy end and insert
        n_count1 = abs(ref2_start - ref1_end)/2
        n_count2 = n_count1
        if abs(ref2_start - ref1_end) % 2 == 1:
            n_count2 += 1
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(n_count2) + "N" + str((cg2_start - cg1_end)) + "I" + str(n_count1) + "N" + cg2_cleaned


# Gets the start and end postions of any any hard and soft clips on a read
# if both at one end combine them together into one segment
def get_read_pos(record, front):
    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)
    if not front:
        cigars = cigars[::-1]
    # Search until we find a softclip. Will happen at start or end. May have a hard clip before
    start = 0
    end = 0
    for cg in cigars:
        if cg.endswith('H'):
            step = int(cg[:cg.find("H")])
            end += step
        elif cg.endswith('S'):
            step = int(cg[:cg.find("S")])
            end += step
        else:
            break
    return end

parser = argparse.ArgumentParser( description='Read the results of mapping reads to a reference genome, output a paf with split mappings combined if needed')
parser.add_argument('--read-to-reference-bam', type=str, required=True)
parser.add_argument('--input-bedpe', type=str, required=False)
parser.add_argument('--output-bedpe', type=str, required=True)
parser.add_argument('--min-insert-size', type=int, default=5000)
parser.add_argument('--reference-window', type=int, default=5000)
parser.add_argument('--min-detected-inclusion-length', type=int, default=50)
parser.add_argument('--reference-gap-minimum', type=int, default=10)
parser.add_argument('--minimum-mapping-qual', type=int, default=20)
args = parser.parse_args()

regions = {}
if args.input_bedpe:
    with open(args.input_bedpe, 'r') as in_bedpe:
        for line in in_bedpe:
            row = line.strip().split('\t')
            regions[row[0]+":"+str(int(row[1])-args.reference_window)+"-"+str(int(row[1])+args.reference_window)] = 1
            regions[row[3]+":"+str(int(row[4])-args.reference_window)+"-"+str(int(row[4])+args.reference_window)] = 1
else:
    # Iterate over the full bam
    header = pysam.AlignmentFile(args.read_to_reference_bam).header
    for sq in header['SQ']:
        length = int(sq['LN'])
        for i in range(1, length, 100000):
            regions[sq['SN']+":"+str(i)+"-"+str(i+100000)] = 1

with open(args.output_bedpe, 'w') as out_merged:
    sam_reader = pysam.AlignmentFile(args.read_to_reference_bam)
    for region in regions:
        print(region)
        tmp_sam_reader = sam_reader.fetch(region=region)
        read_to_reference_alignments = defaultdict(list)
        for record in tmp_sam_reader:
            read_to_reference_alignments[record.query_name].append(record)
        for read_id in read_to_reference_alignments:
            records = sorted(read_to_reference_alignments[read_id], key = lambda x: (x.reference_id, x.reference_start))
            # Assume that records has been sorted by chromsosme and position
            if len(records) > 1:
                # Check for split reads
                record_1 = records[0]
                i = 1
                while i < len(records):
                    record_2 = records[i]
                    if (record_1.query_name != record_2.query_name or 
                        record_1.reference_name != record_2.reference_name or
                        record_1.is_reverse != record_2.is_reverse or
                        (record_1.reference_length + record_2.reference_length) > 1.2* get_read_length(record_1.cigarstring) or
                        (abs(record_1.reference_end - record_2.reference_start) > int(args.reference_gap_minimum) and 
                        abs(record_2.reference_end - record_1.reference_start) > int(args.reference_gap_minimum)) or
                        record_1.mapping_quality < args.minimum_mapping_qual or record_2.mapping_quality < args.minimum_mapping_qual):
                        record_1 = record_2
                        i += 1
                        continue
                    read_length = get_read_length(record_1.cigarstring)
                    if abs(record_1.reference_end - record_2.reference_start) <= int(args.reference_gap_minimum):
                        insert_start = read_length - get_read_pos(record_1, False)
                        insert_end = get_read_pos(record_2, True)
                        if abs(insert_start - insert_end) >= args.min_insert_size:
                            start = min(record_1.reference_end, record_2.reference_start)
                            end = max(record_1.reference_end, record_2.reference_start)
                            out_merged.write(record_1.reference_name+"\t"+str(start)+"\t"+str(end)+"\t"+record_1.query_name+"\t"+str(abs(insert_start - insert_end))+"\n")
                    elif abs(record_2.reference_end - record_1.reference_start) <= int(args.reference_gap_minimum):
                        insert_start = read_length - get_read_pos(record_2, False)
                        insert_end = get_read_pos(record_1, True)
                        if abs(insert_start - insert_end) >= args.min_insert_size:
                            start = min(record_2.reference_end, record_1.reference_start)
                            end = max(record_2.reference_end, record_1.reference_start)
                            out_merged.write(record_1.reference_name+"\t"+str(start)+"\t"+str(end)+"\t"+record_1.query_name+"\t"+str(abs(insert_start - insert_end))+"\n")
                    i += 1
                    if i < len(records):
                        record_1 = records[i]
                        i += 1

            
