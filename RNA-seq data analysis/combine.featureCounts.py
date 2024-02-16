#!usr/bin/python3
#coding:utf-8

import sys,re
if len(sys.argv)<2:
    print("usage: python3 combine.reads.count.py out.file input.files.....")
    sys.exit(-1)
else:
    file_dict={}
    out=open(sys.argv[1],'w')
    keys={}
    files=[]
    for i in sys.argv[2:len(sys.argv)]:
        fh=open(i,'r')
        i=re.split("\/",i)[-1]
        i=re.sub(".featureCounts.*","",i)
        files.append(i)
        for line in iter(fh):
            if not line.startswith("Geneid") and not line.startswith("#"):
                line_info=re.split("\t",line.strip())
                keys[line_info[0]]=1
                file_dict[(i,line_info[0])]=line_info[6]
    out.write("genes")
    for f in files:
        out.write("\t"+f)
    out.write("\n")
    for k in keys.keys():
        out.write(k)
        for f in files:
            if (f,k) in file_dict.keys():
                out.write("\t"+str(file_dict[(f,k)]))
            else:
                out.write("\t0")
        out.write("\n")
