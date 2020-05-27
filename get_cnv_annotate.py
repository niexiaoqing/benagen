#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
f = open(sys.argv[1],'r')
list1 = f.readlines()
list1 = [i.strip() for i in list1]
list1_1 = [i.split()[0] for i in list1]
f.close()
my_set = set(list1_1)

f = open(sys.argv[2],'r')
list2 = f.readlines()
f.close()
list2_1 = [i.split()[0] for i in list2][1:]
my_set1 = set(list2_1)
or_set = my_set - my_set1
and_set = my_set & my_set1


list3 = []
for i in list2:
    if i.split()[0] in and_set:
        list3.append(i)
list3.sort(key=lambda x:x.split()[0])

list4 = []
for i in list1:
    if i.split()[0] in and_set:
        list4.append(i)
list4.sort(key=lambda x:x.split()[0])

list5 = []
for i,j in zip(list3,list4):
    list5.append(j.strip() + "\t" + '\t'.join(i.strip().split("\t")[1:]))

list6 = []
for i in list1:
    if i.split()[0] in or_set:
        list6.append(i + "\t" + "\t".join(["--"]*6)) 

list7 = list6 + list5
list7.sort(key=lambda x:x.split()[1])


print("\t".join(["GENE_ID","REGION",'SOURCE','GO','KEGG','SWISSPROT','PFAM','NR']))
for i in list7:
    print(i)


