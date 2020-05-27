#!/home/benagen/anaconda3/bin/python
import sys
import re
import subprocess
def USAGE():
    if len(sys.argv) < 3:
        print("""usage:
    ch_reads_name.py arg1 arg2
    arg1 is the suffix of reads1
    arg2 is the suffix of reads2
    and you should change dir to rawreads dir as current work dir
""")
        exit(1)

USAGE()


print("aaaaaaaaaaa")

p1 = subprocess.Popen("ls *" + sys.argv[1],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
p1.wait()
p2 = subprocess.Popen("ls *" + sys.argv[2],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
p2.wait()
print(str(p1.stderr.read(),encoding="utf8"))
str1 = str(p1.stdout.read(),encoding="utf8")
print(str(p2.stderr.read(),encoding="utf8"))
str2 = str(p2.stdout.read(),encoding="utf8")
list1 = str1.split() + str2.split()

for i in list1:
    if i.endswith(sys.argv[1])
        idx = i.rindex(sys.argv[1])
        temp = i[0:idx]
        p3 = subprocess.Popen('mv ' + i + ' ' + temp + '_R1.fq.gz',shell=True,stdout=subprocess.PIPE,stderr=SUBPROCESS.pipe)
        p3.wait()
        print(str(p3.stderr.read(),encoding="utf8"))
    if i.endswith(sys.argv[2])
        idx = i.rindex(sys.argv[2])
        temp = i[0:idx]
        p3 = subprocess.Popen('mv ' + i + ' ' + temp + '_R2.fq.gz',shell=True,stdout=subprocess.PIPE,stderr=SUBPROCESS.pipe)
        p3.wait()
        print(str(p3.stderr.read(),encoding="utf8"))













