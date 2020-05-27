#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
import re

def usage():
    print(
"""
第一个参数(sys.argv[1]) : .gff文件，读取
第二个参数(sys.argv[2]) : .cnv文件，读取
第三个参数(sys.argv[2]) : 指定输出文件的路径，写操作
注意：该脚本的代码可能会根据.gff文件的数据结构的差异而修改


"""
)

if len(sys.argv) == 1:
    usage()
    sys.exit()



gff_list = []
f1 = open(sys.argv[1],'r')
for i in f1:
    if (not i.startswith("#")) and i.strip() != '':
        temp = i.strip().split("\t")
        if temp[2] == "gene":
            gff_t = (temp[0],int(temp[3]), int(temp[4]),'_'.join(temp[8].split(';')[0].split('=')[1].split('_')[0:]))
            gff_list.append(gff_t)
f1.close()
gff_list.sort(key=lambda x:x[1])
gff_list.sort(key=lambda x:x[0])

# print(gff_list)
# print('*'*80)

cnvs = []
f2 = open(sys.argv[2],'r')
for i in f2:
    temp = i.strip().split()
    my_region = re.split(':|-',temp[1])
    my_region = (my_region[0] ,int(my_region[1]),int(my_region[2]))
    cnvs.append(my_region)

f2.close()
cnvs.sort(key=lambda x:x[1])   
cnvs.sort(key=lambda x:x[0])   
# print(cnvs)
gff_list1 = [i[0] for i in gff_list]
cnv_list1 = [i[0] for i in cnvs]
my_set1 = set(gff_list1 + cnv_list1)
my_set = sorted(list(my_set1))
# print(my_set)
gff_list = [i for i in gff_list if i[0] in my_set1]
cnvs = [i for i in cnvs if i[0] in my_set1]
a1 = 0
a2 = 0
gff_len = len(gff_list)
cnv_len = len(cnvs)
gene_list = set()

while True:
    if my_set.index(gff_list[a1][0]) > my_set.index(cnvs[a2][0]):
        if a2 == cnv_len - 1:
            break
        a2 += 1 
    elif my_set.index(gff_list[a1][0]) < my_set.index(cnvs[a2][0]):
        if a1 == gff_len - 1 :
            break
        a1 += 1
    else:
        if gff_list[a1][2] < cnvs[a2][1]:
            if a1 == gff_len - 1:
                break
            a1 += 1
        elif gff_list[a1][1] > cnvs[a2][2]:
            if a2 == cnv_len - 1:
                break
            a2+=1
        else:
           # print(cnvs[a2], gff_list[a1])
            gene_list.add(gff_list[a1][3] + '\t' + cnvs[a2][0] + ":" + str(cnvs[a2][1]) + ':' + str(cnvs[a2][2]))
            if a1 == gff_len -1:
                break
            a1 += 1

f3 = open(sys.argv[3],'w')
for i in gene_list:
    f3.write(i + "\n")
f3.close() 


"""        
f3 = open(sys.argv[3],'r',encoding='latin1')
annoate_list = f3.readlines()
f3.close()
for i in annoate_list:
    temp = i.split('\t')
    if temp[0] in gene_list:
        print(i.strip())

"""


