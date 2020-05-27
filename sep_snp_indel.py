#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
import re
pat = re.compile(r'^[ACGT]$')


f = open(sys.argv[1],'r')
f1 = open("all.snp.vcf",'w')
f2 = open('all.indel.vcf','w')



for i in f:
    if i.startswith("#"):
        f1.write(i)
        f2.write(i)
    elif '/2' in i:
        continue
    else:
        temp = i.strip().split()
        if re.search(pat,temp[3]) and re.search(pat,temp[4]):
            f1.write(i)
        else:
            f2.write(i)

f.close()
f1.close()
f2.close()
