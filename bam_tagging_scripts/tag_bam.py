__author__ = "etseng@pacb.com"
"""
Tagging BAM files with phasing info 
"""

import pysam
from csv import DictReader


def main(read_bam, output_bam, namefile, labl, tag_name):
    d = {}
    celltype_info = {}
    for r in DictReader(open(namefile), delimiter=','):
        d[r['id']] = r[labl]

    reader = pysam.AlignmentFile(read_bam, 'rb', check_sq=False)
    f2 = pysam.AlignmentFile(output_bam, 'wb', header=reader.header)
    for r in reader:
        d2 = r.to_dict()
        if r.qname in d:
            d2['tags'].append(tag_name + str(d[r.qname]))
        x = pysam.AlignedSegment.from_dict(d2, r.header)
        f2.write(x)

    f2.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Tagging BAM files with phasing info")
    parser.add_argument("read_bam", help="Aligned BAM file that be tagged")
    parser.add_argument("output_bam", help="Output tagged BAM filename")
    parser.add_argument("namefile", help="[Optional] Comma-delimited namefile info CSV, must have column 'bc' and 'namefile'")
    parser.add_argument("labl", help="[Optional] string for label")
    parser.add_argument("--tag_name", help="[Optional] string to append labels to")

    args = parser.parse_args()
    main(args.read_bam, args.output_bam, args.namefile, args.labl, args.tag_name)
