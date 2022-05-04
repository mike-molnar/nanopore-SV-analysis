import sys
import argparse
import vcf

# Function to return the other half of a translocation call
def get_bnd(alt):
    if alt[0] in [']', '[']:
        chr2 = alt.split(':')[0][1:]
        pos2 = int(alt.split(':')[1][:-2])
    else:
        chr2 = alt.split(':')[0][2:]
        pos2 = int(alt.split(':')[1][:-1])

    return chr2, pos2

parser = argparse.ArgumentParser( description='Convert a structural variant VCF file to a BEDPE file.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

in_fh = vcf.Reader(open(args.input, 'r'))
out_fh = open(args.output, 'w')
out_fh.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_type\tsv_length\tsupport\n")

for record in in_fh:
    try:
        sv_length = abs(record.INFO["SVLEN"])
    except:
        sv_length = 0

    try:
        sv_support = record.INFO["SUPPORT"]
    except:
        sv_support = record.INFO["RE"]
        
    if record.INFO['SVTYPE'] in ["DEL", "INS", "INV", "DUP", "DUP:TANDEM", "DUP:INT"]:
        chr2 = record.CHROM
        start_2 = record.INFO["END"]-1
        end_2 = record.INFO["END"]
    else:
        chr2, pos2 = get_bnd(str(record.ALT[0]))
        start_2 = pos2-1
        end_2 = pos2
        
    out_fh.write("{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{svtype}\t{length}\t{support}\n".format(
            chrom1 = record.CHROM,
            start1 = record.POS-1,
            end1 = record.POS,
            chrom2 = chr2,
            start2 = start_2,
            end2 = end_2,
            svtype = record.INFO['SVTYPE'],
            length = sv_length,
            support = sv_support))
            
out_fh.close()