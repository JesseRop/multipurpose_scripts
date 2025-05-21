import numpy as np
import os
import pysam
import argparse
import errno


def main():
    parser = argparse.ArgumentParser(description="Splits BAM file into two BAM files, putting p reads into one BAM and (1-p) into the other")
    parser.add_argument('-b', '--bam-file', dest='bam_file', required=True, help='Input BAM file')
    parser.add_argument('-o1', '--output1', dest='output1', required=False, help='Output BAM file #1')
    parser.add_argument('-o2', '--output2', dest='output2', required=False, help='Output BAM file #2')
    parser.add_argument('-d', '--output-dir', dest='output_dir', required=False, help='Output directory')
    parser.add_argument('-p', '--parameter', dest='parameter', type=float, required=True, help='Parameter p of binomial distribution')
    args = parser.parse_args()

#    if not os.path.isfile(args.bam_file):
#        print "ERROR: BAM file does not exist!"
#        exit(-1)

#    if args.parameter < 0. or args.parameter > 1.:
#        print "ERROR: parameter must be float value in range 0.0 <= p <= 1.0"
#        exit(-1)

    if not args.output1:
        output1_filename = os.path.splitext(args.bam_file)[0] + "_" + str(args.parameter) + os.path.splitext(args.bam_file)[1]
    else:
        output1_filename = args.output1
    if not args.output2:
        output2_filename = os.path.splitext(args.bam_file)[0] + "_" + str(1.0 - args.parameter) + os.path.splitext(args.bam_file)[1]
    else:
        output2_filename = args.output2

    # if user specified output dir
    if args.output_dir:  
        # if it already exists, we join paths
        if os.path.isdir(args.output_dir):
            output1_filename = os.path.join(args.output_dir, output1_filename)
            output2_filename = os.path.join(args.output_dir, output2_filename)
        # else try and create, then join
        else:
            try:
                os.makedirs(args.output_dir)
                output1_filename = os.path.join(args.output_dir, output1_filename)
                output2_filename = os.path.join(args.output_dir, output2_filename)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
                    
    # open input BAM file for reading
    bam_f = pysam.AlignmentFile(args.bam_file, "rb")
    # open output BAM files for writing 
    output1_f = pysam.AlignmentFile(output1_filename, "wb", template=bam_f)
    output2_f = pysam.AlignmentFile(output2_filename, "wb", template=bam_f)

    num_reads = 0
    num_paired_reads = 0
    num_unpaired_reads = 0
    num_mate_missing = 0
    num_file1 = 0
    num_file2 = 0
    bam_iter = bam_f.fetch(until_eof=True)
    # for each read in BAM file
    for read in bam_iter:
        # sample binomial to choose which output BAM file to write to
        coin_flip = np.random.binomial(1, args.parameter)
        # if paired, get matching pair and write to output
        if read.is_paired and read.is_proper_pair:
            # if this is the first read in the pair, get second read
            if read.is_read1 and not read.mate_is_unmapped:
                bam_iter_tmp = bam_iter
                # get the read's mate
                try:
                    read2 = bam_f.mate(read)
                    if coin_flip:
                        output1_f.write(read)
                        output1_f.write(read2)
                        num_file1 += 2
                    else:
                        output2_f.write(read)
                        output2_f.write(read2)
                        num_file2 += 2
                    num_paired_reads += 2
                    num_reads += 2
                # unless the mate is unavailable for some reason
                except ValueError:
                    if coin_flip:
                        output1_f.write(read)
                        num_file1 += 1
                    else:
                        output2_f.write(read)
                        num_file2 += 1
                    num_mate_missing += 1
                    num_reads += 1
                    continue
                finally:
                    bam_iter = bam_iter_tmp
        # else read is not paired, just throw it into a file randomly
        else:
            if coin_flip:
                output1_f.write(read)
                num_file1 += 1
            else:
                output2_f.write(read)
                num_file2 += 1

            num_unpaired_reads += 1
            num_reads += 1
#       if num_reads % 100000 == 0:
#            print "Number of reads processed:", num_reads
            
    bam_f.close()
    output1_f.close()
    output2_f.close()
#    print "Total reads processed:", num_reads
#    print "Number paired reads:", num_paired_reads, "\t(" + str(100.*num_paired_reads / float(num_reads)) + "%)"
#    print "Number unpaired reads:", num_unpaired_reads, "\t(" + str(100.*num_unpaired_reads / float(num_reads)) + "%)"
#    print "Number paired reads missing their mate:", num_mate_missing, "\t(" + str(100.*num_mate_missing / float(num_reads)) + "%)"
#    print "Number of reads in first file:", num_file1, "\t(" + str(100.*num_file1 / float(num_reads)) + "%)"
#    print "Number of reads in second file:", num_file2, "\t(" + str(100.*num_file2 / float(num_reads)) + "%)"

if __name__ == "__main__":
    main()
