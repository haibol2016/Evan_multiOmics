## Step 0: merge sequence data from the same library

This is only necessary if the same library has been sequenced in
multiple lanes. After this step, each sample/library will have two fastq
files for R1 and R2.

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=2000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-192]%50"
    #BSUB -q short   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs
    midir -p data

    L1_R1=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L1_1.fq.gz" |sort `)
    L4_R1=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L4_1.fq.gz" |sort `)

    L1_R2=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L1_2.fq.gz" |sort `)
    L4_R2=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L4_2.fq.gz" |sort `)

    L1=(${L1_R1[@]} ${L1_R2[@]})
    L4=(${L4_R1[@]} ${L4_R2[@]})

    L1_R1_name=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L1_1.fq.gz" |sort | perl -p -e 's{.+/}{}'`)
    L1_R2_name=(`find /pi/thomas.fazzio-umw/Evan/worm_scCoT/20240612_worm_scCoT/fastq/usftp21.novogene.com/01.RawData -name "*L1_2.fq.gz" |sort | perl -p -e 's{.+/}{}'`)

    name=(${L1_R1_name[@]} ${L1_R2_name[@]})
    cat ${L1[$i]} ${L4[$i]} > data/${name[$i]}

## Step 2: Perform QC on raw reads using FastQC

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=2000] # ask for memory 5G
    #BSUB -W 2:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-192]%50"
    #BSUB -q short   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    module load fastqc/0.11.9
    in=data

    fastq=(`ls $in/*.gz`)
    out=results/001.fastqc.out
    mkdir -p $out

    fastqc -o $out  -t 4 --extract ${fastq[$i]}

## Step 3: Generate an integrative QC report using MultiQC

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=4000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1]"
    #BSUB -q short   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(89472)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate multiqc1.14


    in_dir=(results/001.fastqc.out) 
    out_dir=results/002.MultiQC.out
    project=(Raw.data)

    mkdir -p ${out_dir}

    ## using virtual environment

    multiqc --filename  ${project[$i]}.multiQC  --outdir  ${out_dir}   ${in_dir[$i]}

## Step 4: Extract cell barcodes using UMItools

See UMItools document: [UMItools
extract](https://umi-tools.readthedocs.io/en/latest/reference/extract.html).
Here the “regex” method is used to extract cell barcode, with one cell
barcode treated as UMI to make **UMItool extract** command runnable.

For specifying the cell barcode patterns, please refer to documentation
of the python package **regex**: [regex fuzzy
match](https://github.com/mrabarnett/mrab-regex?tab=readme-ov-file#approximate-fuzzy-matching-hg-issue-12-hg-issue-41-hg-issue-109).

A fuzzy regex specifies which types of errors are permitted, and,
optionally, either the minimum and maximum or only the maximum permitted
number of each type. (You cannot specify only a minimum.)

The 3 types of error are:

-   Insertion, indicated by “i”
-   Deletion, indicated by “d”
-   Substitution, indicated by “s” In addition, “e” indicates any type
    of error.

The fuzziness of a regex item is specified between “{” and “}” after the
item. The fuzziness should be determined based on the average base
quality per read. After extraction, the cell barcodes and so-called UMI,
which is Antibody barcode here, will be appended to the first part of
the names of R1 and R2 in the format “\_CCAGTTCAAGATGT\_TCTCCGGA”, ie,
“\_CellBarcode1CellBarcode2\_UMI”, as shown follows:

@A00742:875:HHCNLDSXC:1:1101:18186:1047\_CCAGTTCAAGATGT\_TCTCCGGA
1:N:0:TAGCTT+TGGTCA
CAGCTGGCCCTGGAGGTCCTGGTGGTCCCTGACCTCCACCGTGTGCTGGGGCT +
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFF

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=2000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-96]%50"
    #BSUB -q long   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs


    module load umi_tools/1.1.2
    r1=(data/*_1.fq.gz)
    name=(`ls data/*_1.fq.gz | perl -p -e 's{.+/(.+?)_1.fq.gz}{$1}'`)
    r2=(data/*_2.fq.gz)
    out=data.cellbarcode/ab_barcode_demultiplex
    mkdir -p $out

    #ref:https://github.com/mrabarnett/mrab-regex

    umi_tools extract --extract-method=regex \
        --bc-pattern="(?P<cell_1>^.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){s<=4,i<=2,d<=2}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCAAAGAAAGATGTGTATAAGAGACAG){e<=8}.*" \
        --bc-pattern2="(?P<discard_3>GTCTCGTGGGCTCGGCTGTCCCTGTCC){s<=3,i<=2,d<=2}(?P<umi_1>TCTCCGGA|AATGAGCG|GGAATCTC){s<=1}(?P<discard_4>AGATGTGTATAAGAGACAG){e<=6}.*"  \
        --read2-in=${r2[$i]} --filtered-out $out/${name[$i]}_UNKNOWN_1.fq.gz \
        --filtered-out2  $out/${name[$i]}_UNKNOWN_2.fq.gz \
        --read2-out=$out/${name[$i]}_2.fq.gz -L logs/${name[$i]}.extract.log  \
        -I ${r1[$i]} -S $out/${name[$i]}_1.fq.gz \
        --compresslevel=9

## Step 5: Split fastq file of each sample by antibody-specific barcode using

custom python scripts: “005splitbyab.py”

At this step, the last round cell barcodes–96-well-specific sample
indices–is appended to the concatenated previous rounds’ cell barcodes,
and UMI (antibody- specific barcode) are used to split the fastq file
into antibody-specific fastq files.

Running environment: Any conda environment with Python 3, gzip, and
regex python package installed.

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

Run the python scripts as follows.

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=2000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-96]%50"
    #BSUB -q long   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate cleanuprnaseq    

    in=data.cellbarcode/ab_barcode_demultiplex
    r1=($in/*L1_1.fq.gz)
    name=(`ls data/*_L1_1.fq.gz | perl -p -e 's{.+/(.+?)_L1_1.fq.gz}{$1}'`)
    r2=($in/*L1_2.fq.gz)
    out=data.ab.split.out
    mkdir -p $out


    python scripts/003.splitbyab.py ${r1[$i]} ${r2[$i]} ${name[$i]} $out

## Step 6: Merge fastq files by antibody

Step 5 generated 3 pairs of fastq files for each input pair of fastq
(96). To reduce file numbers, we need merge them by antibody.

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=4000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1]"
    #BSUB -q short   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(89472)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    in=data.ab.split.out
    h3k27_1=(`ls $in/*H3K27me3_1.fq.gz | sort`)
    h3k27_2=(`ls $in/*H3K27me3_2.fq.gz |sort`)
    h3k4_1=(`ls $in/*H3K4me3_1.fq.gz | sort`)
    h3k4_2=(`ls $in/*H3K4me3_2.fq.gz  | sort`)
    pol_1=(`ls $in/*Pol2-S2P_1.fq.gz  | sort`)
    pol_2=(`ls $in/*Pol2-S2P_2.fq.gz  | sort`)

    out=data/00.antibody.split.out
    mkdir -p ${out}

    cat ${h3k27_1[@]} > $out/H3K27me3_S1_L001_R1_001.fq.gz
    cat ${h3k27_2[@]} >  $out/H3K27me3_S1_L001_R2_001.fq.gz
    cat ${h3k4_1[@]} >  $out/H3K4me3_S1_L001_R1_001.fq.gz
    cat ${h3k4_2[@]} >  $out/H3K4me3_S1_L001_R2_001.fq.gz
    cat ${pol_1[@]} >   $out/PolII_S1_L001_R1_001.fq.gz
    cat ${pol_2[@]} >  $out/polII_S1_L001_R2_001.fq.gz

## Step 7: Index **C. elegans** reference genome using BWA

    #!/bin/bash

    #BSUB -n 1  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=8000] # ask for memory 8G per core
    #BSUB -W 72:00 #limit the job to be finished in 4 hours
    #BSUB -J "bwa_index"
    #BSUB -q long # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"

    mkdir -p logs

    module load bwa/0.7.17

    gunzip -k docs/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

    fasta=docs/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa

    bwa index -p docs/WBcel235  ${fasta}

## Step 8: Map reads to the reference genome using BWA mem

    #!/bin/bash

    #BSUB -n 8  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=8000] # ask for memory 5G
    #BSUB -W 8:00 #limit the job to be finished in 12 hours
    #BSUB -J "bwa[1-3]"
    #BSUB -q long  # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(9390569)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    module load bwa/0.7.17 
    module load samtools/1.16.1 


    in=data/00.antibody.split.out/

    R1=(`ls $in/*R1_001.fq.gz`)
    R2=(`ls $in/*R2_001.fq.gz`)


    name=(`ls  $in/*_R1_001.fq.gz | perl -p -e 's{.+/(.+?_.+?)_.+?R1.+?.gz}{$1}g'`)
    pu=(`ls  $in/*_R1_001.fq.gz | perl -p -e 's{.+/.+?_.+?_(.+?)_R1.+?.gz}{$1}g'`)
    out=results/003.2.bwa.out/
    mkdir -p ${out}

    fasta=docs/WBcel235

    time bwa mem -M -t 8  \
      -R "@RG\tID:${name[$i]}\tPL:ILLUMINA\tSM:${name[$i]}\tPU:${pu[$i]}\tLB:${name[$i]}"  $fasta  ${R1[$i]} ${R2[$i]} | \
         samtools sort  -l 5 -m 8G -o $out/${name[$i]}_${pu[$i]}.sort.bam -O BAM -@ 8  -

    samtools index $out/${name[$i]}_${pu[$i]}.sort.bam

## Step 9: Generate fragments file using Sinto

See documentation for **sinto fragments**: [Sinto
fragments](https://timoast.github.io/sinto/basic_usage.html#create-scatac-seq-fragments-file).

    #!/bin/bash

    #BSUB -n 4  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=8000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-3]%50"
    #BSUB -q long   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs
    module load sinto

    in=results/003.2.bwa.out
    bam=(`ls $in/*.bam`)
    name=(`ls $in/*.bam | perl -p -e 's{.+/(.+?).sort.bam}{$1}'`)
    out=results/004.fragments.out

    mkdir -p $out

    sinto fragments -b ${bam[$i]} -f $out/${name[$i]}.fragments -m 20 -p 4 \
                       --barcode_regex "[ATCG]+$"  --use_chrom "^[IVX]" \
                       --max_distance 2000 --min_distance 10 \
                       --chunksize 1000000 --shift_plus 4 \
                       --shift_minus -5 --collapse_within

## Step 10: Sort and index fragments file using samtools

For fragments file format, see [10x
Genomic](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/fragments-file).

For library quality check up, we can visualize the fragments files in
[IGV](https://igv.org/app/) before go to Step 11

    #!/bin/bash

    #BSUB -n 4  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=8000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-3]%50"
    #BSUB -q long   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs
    module load samtools


    in=results/004.fragments.out
    frag=(`ls $in/*.fragments`)
    name=(`ls $in/*.fragments | perl -p -e 's{.+/(.+?).fragments}{$1}'`)
    out=results/004.fragments.out

    mkdir -p $out

    sort -k 1,1 -k2,2n  ${frag[$i]} |bgzip -c -f > $out/${name[$i]}.fragments.tsv.gz
    tabix -p bed  $out/${name[$i]}.fragments.tsv.gz 

## Step 11: Cell clustering, peak calling using ArchR

Here we treat cut&tag data as scATAC-seq data, and perform analysis
using [ArchR](https://www.archrproject.com/). However, ArchR by default
only support analysis for human and mouse data. To use ArchR for any
species scATAC-seq data, we developed and R package,
[GenomePal](https://github.com/haibol2016/GenomePal).

It is not easy to install ArchR, so we build a docker image available
from {Docker Hub\](hukai916/scatacpipe\_downstream:0.2.1), which can be
pulled and converted into a singularity image and run interactively.

I will provide further details if you think you can run this step.

## Step 12: Split fragments file by Cluster using custom Python scripts

    !/usr/bin/env python

    """
    Split bed files into clusters given cluster file and fragment file.

    Usage:
    python split_bed.py Cluster_xxx.tsv Fragment.bed.gz

    Cluster_xxx.tsv file must consist of two columns (first column is cell barcode, second is the assigned cluster), separated by tab.
    Fragment.bed.gz must be in .bed or .bed.gz format, and the last second column must be cell barcode.
    """

    import sys
    import csv
    import gzip
    import os
    import shutil
    from pathlib import Path

    cluster_file    = sys.argv[1]
    fragment_file   = sys.argv[2]
    out_dir          = sys.argv[3] 
    cluster_dict    = {}

    with open(cluster_file) as tsv_file:
        tsv = csv.reader(tsv_file, delimiter='\t')
        next(tsv)
        for row in tsv:
            cluster_dict[row[-1]] = row[-2]

    clusters = set(x for x in cluster_dict.values())

    # create output folder using cluster_file name:
    out_dir = os.path.join(out_dir, "split_" + Path(cluster_file).stem)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)


    def split_bed(f, cluster_dict, output_cluster_dict):
        for i, line in enumerate(f):
            if i % 1000000 == 0:
                print(str(int(i / 1000000) * 1000000) + " reads processed ...")
            if not line.startswith("#"): # in case some software generated bed file starts with comment lines
                try: # some barcode are not in the cluster_dict
                    barcode = line.split()[-2]
                    if cluster_dict[barcode] in output_cluster_dict:
                        output_cluster_dict[cluster_dict[barcode]].append(line)
                    else:
                        output_cluster_dict[cluster_dict[barcode]] = [line]
                except:
                    continue

    def save_file(out_dir, output_cluster_dict):
        print("Saving to output ...")
        for cluster in output_cluster_dict:
            output_file = out_dir + "/cluster_" + cluster + ".txt"
            with open(output_file, "w") as res:
                res.write("".join(output_cluster_dict[cluster]))

    # loop over fragment file and output to corresponding cluster file.
    output_cluster_dict = {}

    if fragment_file.endswith(".gz"):
        with gzip.open(fragment_file, mode='rt') as f:
            split_bed(f, cluster_dict, output_cluster_dict)
        save_file(out_dir, output_cluster_dict)
    elif fragment_file.endswith(".bed"):
        with open(fragment_file, mode='rt') as f:
            split_bed(f, cluster_dict, output_cluster_dict)
        save_file(out_dir, output_cluster_dict)
    else:
        print("Pls supply either .bed or .bed.gz file.")
        exit()

    print("Done!")

    #!/bin/bash

    #BSUB -n 4  # minmal numbers of processors required for a parallel job
    #BSUB -R rusage[mem=8000] # ask for memory 5G
    #BSUB -W 4:00 #limit the job to be finished in 12 hours
    #BSUB -J "fastQC[1-2]%50"
    #BSUB -q long   # which queue we want to run in
    #BSUB -o logs/out.%J.%I.txt # log
    #BSUB -e logs/err.%J.%I.txt # error
    #BSUB -R "span[hosts=1]" # All hosts on the same chassis"
    ##BSUB -w "done(5423513)"

    i=$(($LSB_JOBINDEX- 1))
    mkdir -p logs

    source /home/haibo.liu-umw/miniconda3/etc/profile.d/conda.sh
    conda activate cleanuprnaseq

    in=results/004.fragments.out
    frag=(results/004.fragments.out/{H3K27me3_S1_L001.fragments.bed,H3K4me3_S1_L001.fragments.bed})
    out=results/005.split.bed

    mkdir -p $out

    cluster=(docs/H3K27me3.cluster.info.txt docs/H3K4me3.cluster.info.txt)


    python scripts/010.split.bed.py ${cluster[$i]}  ${frag[$i]}  $out
