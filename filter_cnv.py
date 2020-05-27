#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
import subprocess



p1 = subprocess.Popen("ls *.cnv",shell=True,stdout=subprocess.PIPE)
p1.wait()
s1 = p1.stdout.read()
s2 = str(s1,encoding="utf8")

if '.set.cnv' in s2:
    print("*.set.cnv already exists")
    exit(1)

list1 = str(s1,encoding="utf8").strip().split()

for i in list1:
    temp = i.strip().split('.')
    temp1 = '.'.join(temp[0:-1])
    f1 = open(temp1 + '.set.cnv','w')
    f1.write("Variation_type\tPOS\tCNV_length\te-value\n")
    f = open(temp1 + '.cnv','r')
    for j in f:
        temp2 = j.strip().split()
        f1.write("\t".join(temp2[0:3]) + "\t" + temp2[4] + "\n")
    f.close()
    f1.close()















