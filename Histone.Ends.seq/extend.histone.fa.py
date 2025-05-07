#!usr/bin/python3
#coding:utf-8
import re,sys
if len(sys.argv)<2:
    print("python3 extend.histone.fa.py histone.gene.list hg38.fa histone.gtf pc_transcript.fa histone.extended.fa extended_length histone.gene.loc.bed")
    sys.exit(-1)
else:
    histone_gene=open(sys.argv[1],'r')
    from Bio import SeqIO
    import subprocess as sp
    histone_gtf=open(sys.argv[3],'r')
    out=open(sys.argv[5],'w')
    ext=int(sys.argv[6])
    out2=open(sys.argv[7],'w')
    histone_gene_dict,transcript,transcript_loc={},{},{}
    for line in iter(histone_gene):
        line_info=re.split("\t",line.strip())
        histone_gene_dict[line_info[7]]=1               # gene ID => 1
    for line in iter(histone_gtf):
        line_info=re.split("\t",line.strip())
        if line_info[2] == "transcript":
            transcript_id=re.split("\"",line_info[8])[3]
            transcript_loc[transcript_id]=[line_info[0],line_info[3],line_info[4],line_info[6]]
    records=list(SeqIO.parse(sys.argv[4],"fasta"))
    longest_transcript_CDS_name,transcript,longest,kept_transcript,longest_transcript_CDS={},{},{},{},{}
    for i in records:
        seq=str(i.seq)
        gene_id=re.split("\|",i.id)[1]
        if gene_id in histone_gene_dict.keys():
            trancript_id=re.split("\|",i.id)[0]
            gene_name=re.split("\|",i.id)[5]
            length=int(re.split("\|",i.id)[6])
            cds_loc=re.findall("CDS:\d+-\d+",i.id)[0].replace("CDS:","")
            transcript[trancript_id]=seq
            if gene_id in longest.keys():
                if length >longest[gene_id]:
                    longest[gene_id]=length
                    kept_transcript[gene_id]=trancript_id
                    longest_transcript_CDS[gene_id]=cds_loc
                    longest_transcript_CDS_name[gene_id]=gene_name
            else:
                longest[gene_id]=length
                kept_transcript[gene_id]=trancript_id
                longest_transcript_CDS[gene_id]=cds_loc
                longest_transcript_CDS_name[gene_id]=gene_name
    for k in histone_gene_dict.keys():
        transcript_seq=transcript[kept_transcript[k]]
        transcript_location=transcript_loc[kept_transcript[k]]
        out2.write(transcript_location[0]+"\t"+str(transcript_location[1])+"\t"+str(transcript_location[2])+"\t"+transcript_location[3]+"\t"+longest_transcript_CDS_name[k]+"\t"+k+"\t"+kept_transcript[k]+"\n")
        exteded_seq=""
        if transcript_location[3] == "+":
            extended=int(transcript_location[2])+ext
            region=transcript_location[0]+":"+str(int(transcript_location[2])+1)+"-"+str(extended)
            if ext >0:
                exteded_seq=re.split("\n",sp.check_output(["samtools","faidx","hg38.fa",region]).decode("utf-8"),1)[1].replace("\n","")
        else:
            extended=int(transcript_location[1])-ext
            region=transcript_location[0]+":"+str(extended)+"-"+str(int(transcript_location[1])-1)
            if ext >0:
                exteded_seq=re.split("\n",sp.check_output(["samtools","faidx","hg38.fa",region,"-i"]).decode("utf-8"),1)[1].replace("\n","")
        out.write(">"+longest_transcript_CDS_name[k]+"\n"+transcript[kept_transcript[k]]+exteded_seq+"\n")
