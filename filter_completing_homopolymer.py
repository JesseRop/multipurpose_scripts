import argparse
import vcf
import pyfaidx

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

fasta = pyfaidx.Fasta(args.fasta)
vcf_reader = vcf.Reader(filename=args.vcf)

vcf_writer = vcf.Writer(open(args.output,'w'), template=vcf_reader)

for record in vcf_reader:
    pos = record.POS - 1
    chrom = record.CHROM
    is_homopolymer = True
    hp_base = None
    hp_length = 2
    print("new rec")
    print(record)
    print(fasta[chrom][pos-hp_length:pos+hp_length+1])
    print(fasta[chrom][pos],record.REF)
    print("reverse")
    if record.is_indel:
        is_homopolymer = False # too complicated to tell
    for index in range(pos-hp_length,pos):
        print(hp_base, fasta[chrom][index])
        if hp_base == None:
            hp_base = fasta[chrom][index]
        elif not(hp_base == fasta[chrom][index]):
            is_homopolymer = False
            break
    print("forward")
    for index in range(pos+1, pos+hp_length+1):
        print(hp_base, fasta[chrom][index])
        if not(hp_base == fasta[chrom][index]):
            is_homopolymer = False
            break
    if is_homopolymer:
        ref_base = fasta[chrom][pos]
        if not(ref_base == hp_base):
            if str(record.ALT[0]) == hp_base: 
                record.FILTER = "hp_completion"
    print(is_homopolymer, record.FILTER)
    vcf_writer.write_record(record)

