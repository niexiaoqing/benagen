#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
f1 = open(sys.argv[1],'r')
list1 = f1.readlines()
f1.close()


list1 = [i.strip() for i in list1]
list1[0] = list1[0] + "\tORG\tGO\tKEGG\tSwissprot\tPfam\tNR"




f2 = open(sys.argv[2],'r')
list2 = f2.readlines()
f2.close()



list2 = [i.strip() for i in list2]
list2 = list2[1:]
list3 = [i.split()[0] for i in list2]
my_set = set(list3)


len1 = len(list1) 
for i in range(1,len1):
    temp = list1[i].split()[5]
    if temp in my_set:
        xx = list2[list3.index(temp)]
        xx1 = "\t".join(xx.split("\t")[1:])
        list1[i] = list1[i] + "\t" + xx1
    else:
        xx = "\t--\t--\t--\t--\t--\t--"
        list1[1] = list1[i] + xx

f2 = open(sys.argv[3],'w')
for i in list1:
    f2.write(i + "\n")

f2.close()

    














