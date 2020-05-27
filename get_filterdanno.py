#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
f = open("all.snpeff.txt",'r')
f1 = open("all.HIGH.snpeff",'w')
f2 = open("all.LOW.snpeff",'w')
f1.write("CHROM\tPOS\tREF\tALT\tMutation_position_type\tgeneid\tprotein\n")
f2.write("CHROM\tPOS\tREF\tALT\tMutation_position_type\tgeneid\tprotein\n")




for i in f:
    if i.startswith("#"):
        continue
    else:
        temp = i.strip().split()
        temp1 = temp[7].split('|')
        if temp1[2] == 'HIGH':
            print(temp[0],temp[1],temp[3],temp[4],temp1[1],temp1[4],temp1[10],sep="\t",file=f1)
        if temp1[2] == 'LOW':
            print(temp[0],temp[1],temp[3],temp[4],temp1[1],temp1[4],temp1[10],sep="\t",file=f2)
f1.close()
f2.close()
f.close()


