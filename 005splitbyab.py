import gzip
import sys
import regex
import os

read1=sys.argv[1]
read2=sys.argv[2]
sample_name = sys.argv[3]
out_dir = sys.argv[4]

cellbarcode_p=regex.compile(r'(.+?_.{14})_(.{8})(\s+[12]:N:0:(.+?)\+(.+))')

h3k27me3_1_f = gzip.open(os.path.join(out_dir, sample_name + "H3K27me3_1.fq.gz"), "at")
h3k27me3_2_f = gzip.open(os.path.join(out_dir, sample_name + "H3K27me3_2.fq.gz"), "at")

h3k4me3_1_f = gzip.open(os.path.join(out_dir, sample_name + "H3K4me3_1.fq.gz"), "at")
h3k4me3_2_f = gzip.open(os.path.join(out_dir, sample_name +"H3K4me3_2.fq.gz"),"at")

pol2_1_f = gzip.open(os.path.join(out_dir, sample_name + "Pol2-S2P_1.fq.gz"), "at")
pol2_2_f = gzip.open(os.path.join(out_dir, sample_name + "Pol2-S2P_2.fq.gz"),"at")

def read_one_read_from_fastq(file_handle):
    name = file_handle.readline().strip()
    seq = file_handle.readline()
    file_handle.readline()
    quality = file_handle.readline()
    return name,seq,quality

with gzip.open(read1,'rt') as r1_f, gzip.open(read2, 'rt') as r2_f:
    # read 1 first read
    (name_1, seq_1, quality_1) = read_one_read_from_fastq(r1_f)
    # read 1 first read
    (name_2, seq_2, quality_2) = read_one_read_from_fastq(r2_f)
    while name_1 and name_2:
        read1_head_info= cellbarcode_p.search(name_1)
        read2_head_info = cellbarcode_p.search(name_2)
        if read1_head_info and read2_head_info:
            name_1 = "".join(read1_head_info.group(1, 4, 5, 3)) + "\n"
            name_2 = "".join(read2_head_info.group(1, 4, 5, 3)) + "\n"
            ab_barcode = regex.match(r'((TCTCCGGA)|(AATGAGCG)|(GGAATCTC)){s<=2}',
                                     read1_head_info.group(2)).group()
            if ab_barcode == "TCTCCGGA":
                h3k27me3_1_f.writelines([name_1, seq_1, "+\n", quality_1])
                h3k27me3_2_f.writelines([name_2, seq_2, "+\n", quality_2])
            elif ab_barcode == "AATGAGCG":
                h3k4me3_1_f.writelines([name_1, seq_1, "+\n", quality_1])
                h3k4me3_2_f.writelines([name_2, seq_2, "+\n", quality_2])
            elif ab_barcode == "GGAATCTC":
                pol2_1_f.writelines([name_1, seq_1, "+\n", quality_1])
                pol2_2_f.writelines([name_2, seq_2, "+\n", quality_2])
            # read 1 first read
            (name_1, seq_1, quality_1) = read_one_read_from_fastq(r1_f)
            # read 1 first read
            (name_2, seq_2, quality_2) = read_one_read_from_fastq(r2_f)

# close file handles
h3k27me3_1_f.close()
h3k27me3_2_f.close()
h3k4me3_1_f.close()
h3k4me3_2_f.close()
pol2_1_f.close()
pol2_2_f.close()
