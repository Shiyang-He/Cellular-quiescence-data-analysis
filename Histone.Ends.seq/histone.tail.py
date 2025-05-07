#!usr/bin/python3
#coding:utf-8
import sys,re,pysam
if len(sys.argv)<2:
    print("usage: python3 filter.bam.file.py sort.by.name.bam dedup.bam tails.table\n")
    sys.exit(-1)
else:
    samfile = pysam.AlignmentFile(sys.argv[1],'rb')
    out_sam=pysam.AlignmentFile(sys.argv[2],'wb',template=samfile)
    tail_table=open(sys.argv[3],'w')
    reads_to_keep,umis={},{}
    for read in samfile.fetch(until_eof=True):
        if read.is_paired and read.is_mapped:
            if re.search("#1:N:0:",read.query_name):
                read_1_id,umi=re.split("_",read.query_name)
                umi_key=(umi,read.reference_name,read.reference_start)
                if umi_key not in umis.keys():
                    out_sam.write(read)
                    umis[umi_key]=1
                    reads_to_keep[read_1_id]=1
            elif re.search("#2:N:0:",read.query_name):
                read_2_id=read.query_name.replace("#2:N:0:","#1:N:0:")
                if read_2_id in reads_to_keep.keys():
                    out_sam.write(read)
                    tail,splice="","" 
                    if re.search("\d+S$",read.cigarstring):## the reads has tails
                        splice=int(re.search("\d+S$",read.cigarstring)[0].replace("S",""))
                        tail=read.query_sequence[-splice:]
                    elif re.search("\d+M$",read.cigarstring):## the reads has no tails
                        splice=0
                        tail="*"
                    if splice !="":
                        output_line="\t".join(map(str,[read.query_name,read.query_sequence,tail,splice,read.reference_end+splice,read.reference_name]))
                        tail_table.write(output_line+"\n")
    print(len(reads_to_keep.keys()))

