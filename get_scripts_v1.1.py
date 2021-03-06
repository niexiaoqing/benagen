#!/data/user/niexiaoqing/anaconda3/bin/python3
import os
import subprocess
import re
import time
pat1 = re.compile(r'#+')

ff = open('./inputs','r')
list11 = ff.readlines()
list12 = []
for i in list11:
    ii = i.strip()
    if ii == '':
        continue
    iii =  re.split(pat1,ii)
    if iii[0] !='':
        iii[0] = iii[0].strip()
        list12.append(iii[0])
ff.close()
list11 = list12
list12=''

dic1 = {}
for i in list11: # GENOME,FASTQ,GFF
    temp = i.split('=')
    dic1[temp[0]] = temp[1]  

R1 = dic1['R_suffix1']
R2 = dic1['R_suffix2']
C1 = dic1['C_suffix1']
C2 = dic1['C_suffix2']
separa = dic1['separation']



GENOME = dic1['GENOME'].split('/')[-1]
pro_path = os.getcwd()
RPATH = pro_path + '/01.rawdata/'
CPATH = pro_path + '/02.cleandata/'
GPATH = pro_path + '/genome/'
LOG = pro_path + '/log/'
path_03 = pro_path + '/03.mapping/'
path_04 = pro_path + '/04.snp_indel/'
path_05 = pro_path + '/05.snp/'
path_06 = pro_path + '/06.indel/'
path_07 = pro_path + '/07.cnv/'
path_08 = pro_path + '/08.annotation/'

genome_path = pro_path + '/data/genomes/'
gff_path = pro_path + '/data/' + '.'.join(GENOME.split('.')[0:-1])
REF_DIR = GPATH + 'Ref_dir/'
BIN = pro_path + '/bin/'


for i in [RPATH,CPATH,GPATH,LOG,path_03,path_04,path_05,path_06,genome_path,gff_path,REF_DIR,BIN]:
    pp1 = subprocess.Popen('[[ -e ' + i + ' ]] || mkdir -p ' + i,shell=True)    
    pp1.wait()
p1 = subprocess.Popen('ls ' + dic1['FASTQ'] + "*gz",shell=True,stdout=subprocess.PIPE)
p1.wait()
ss1 = p1.stdout.read()
print(ss1)
list12 = str(ss1,encoding='utf8').strip().split()
list12 = [i.split('/')[-1] for i in list12]
print(list12)              #######获取文件名basename

for i in list12:
    pp1 = subprocess.Popen('cd ' + RPATH + ' && ln -s ' + dic1["FASTQ"] + i + ' ./',shell=True)
    pp1.wait()
p2 = subprocess.Popen('cd ' + genome_path + '&& ln -s ' + dic1['GENOME'] + ' ' + '.'.join(GENOME.split('.')[0:-1]) + '.fa',shell=True,stderr=subprocess.PIPE)
p2.wait()
print(str(p2.stderr.read(),encoding="utf-8"))
p3 = subprocess.Popen('cd ' + GPATH + '&& ln -s ' + dic1['GENOME'] + ' ./',shell=True,stderr=subprocess.PIPE)
p3.wait()
print(str(p3.stderr.read(),encoding="utf-8"))
p4 = subprocess.Popen('cd ' + gff_path + '&& ln -s ' + dic1['GFF'] + ' genes.gff',shell=True,stderr=subprocess.PIPE)
p4.wait()
print(str(p4.stderr.read(),encoding="utf-8"))
#*************************************** get all dir path and create link file
ggpath = subprocess.Popen('ls ' + GPATH + GENOME, shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
ggpath.wait()

ggpath1 = str(ggpath.stdout.read(),encoding='utf8').strip()
rref_path = '/'.join(GPATH.split('/')[0:-1]) + '/Ref_dir/'
f1 = ''
time.sleep(1)
f = open(ggpath1,'r')
for i in f:
    if i.startswith(">"):
        temp = i.strip().split()[0][1:]
        if f1 != '':
            f1.close()
        f1 = open(rref_path + temp + '.fa','w')
    f1.write(i) 
f.close()
f1.close()
##*****************************************get fastq files in Ref_dir
child1 = subprocess.Popen("find " + RPATH + ' -name "*fq.gz"',shell=True,stdout=subprocess.PIPE)
child1.wait()
str1 = child1.stdout.read()
str1 = str(str1,encoding='utf8').strip()
list1 = str1.split()   # full path of raw_reads(include file name)
list2 = [i.split('/')[-1].split('_')[0] for i in list1]
list2 = []
for i in list1:
    temp = i.split('/')[-1]
    temp1 = temp.split(separa)[0:-1]
    temp2 = separa.join(temp1)
    list2.append(temp2)
list3 = sorted(list(set(list2)))    # sample names(No suffix)

#**************************************
""" 以上代码从fastq文件夹生成包含样本名的文件sample_list.txt             """

script_header = open('software_path.txt','r')
script_header1 = script_header.read()
script_header.close()

#以上代码获得cnvnator,bwa,gatk在运行过程中所需的索引文件。
f = open('./bin/02_run_fastqc.sh','w')
f.write(script_header1 + "\n")
f.write('$fastqc -o ' + RPATH + ' -t 15')
for i in list3:
    f.write(" " + RPATH + i + R1 + ' ' + RPATH + i +R2 )
f.write('\n')
f.write('$fastqc -o ' + CPATH + ' -t 15')
for i in list3:
    f.write(" " + CPATH + i + C1 + ' ' + CPATH + i +C2 )
f.close()
# 获取运行fastqc的shell脚本

f = open('./bin/01_filter.sh','w')
f.write(script_header1 + '\n')
for i in list3:
    s1 = "$soapnuke filter -l 20 -q 0.2 -f $adapter1 -r $adapter2 -1 {0} -2 {1}  -C {2} -D {3} -o {4} \n"
    f.write(s1.format(RPATH + i + R1,RPATH + i + R2,i + C1,i + C2,CPATH[0:-1]))
f.close()
# 数据过滤获得clean_reads
f = open('./bin/03_get_index.sh','w')
f.write(script_header1 + '\n')
f.write("$bwa index " + GPATH  + GENOME + '\n')
f.write("$samtools faidx " + GPATH + GENOME+ '\n')
f.write('java -Xmx2g -jar $create_dict R=' + GPATH + GENOME + ' O=' + GPATH + '.'.join(GENOME.split('.')[0:-1]) + '.dict\n') 
f.close()
# 获取索引文件，位于genome文件夹下
f = open('./bin/04_run_bwa.sh','w')
f.write(script_header1 + '\n')
for i in list3:
    RQ = '"@RG\\tID:{0}\\tLB:{0}\\tPL:illumina\\tPU:{0}\\tSM:{0}"'.format(i)
    f.write('$bwa mem -t 15 -R ' + RQ + ' ' + GPATH +GENOME+' ' +CPATH+i+C1 + ' ' + CPATH + i + C2 + ' -o ' + path_03 + i +'.sam\n') 
f.close()
# **获取运行bwa的shell脚本。
f = open('./bin/05_run_samtools.sh','w')
f.write(script_header1 + '\n')
for i in list3:
    f.write('$samtools sort -@ 10 -O bam  -T ' + i + ' -o ' + path_03 + i + '.sort.bam ' + path_03 + i + '.sam 1>' + LOG + i + '.sam.log 2>' + LOG + i + '.sam.log2\n')
for i in list3:
    f.write('$samtools rmdup -s ' + path_03  + i + '.sort.bam ' + path_03 + i + '.rmdup.bam 1>' + LOG + i + '.rmdup.bam.log 2>' + LOG + i + '.rmdup.bam.log2\n')
for i in list3:
    f.write('$samtools index ' + path_03 + i + '.rmdup.bam\n')
    f.write('$samtools flagstat ' + path_03 + i + '.rmdup.bam >' + LOG + i + '.map_ratio.log\n')
f.close()
# 获取运行samtools的脚本
f = open('./bin/06_run_gatk.sh','w')
f.write(script_header1 + '\n')
f.write('java -jar $gatk3 -T RealignerTargetCreator -nt 75 -R ' +  GPATH + GENOME)
for i in list3:
    f.write(' -I ' + path_03 + i + '.rmdup.bam ')
f.write('-o ' + path_04 + 'all.intervals\n')
f.write('java -jar $gatk3 -T IndelRealigner -R ' + GPATH + GENOME + ' --targetIntervals ' + path_04  + 'all.intervals')
for i in list3:
    f.write(' -I ' + path_03 + i + '.rmdup.bam')
f.write(' -o ' + path_04 + 'all.bam\n')
f.write('$gatk HaplotypeCaller -R ' + GPATH + GENOME + ' -I ' + path_04 + 'all.bam -O '+ path_04 + 'all.gatk.vcf\n')
f.close()
#***********获得运行GATK的脚本
f = open('./bin/07_run_snpeff.sh','w')
f.write(script_header1 + '\n')
f.write('sed \'s#\(data.dir = \)\(.*\)#\\1' + pro_path + '/data/#\' ' + '/data/user/niexiaoqing/tools/snpEff/MODEL>' + pro_path + '/SNPEFF.config\n')
ss = "#{0} genome,version {0}\\n{0}.genome:{0}".format('.'.join(GENOME.split('.')[0:-1]))
f.write('sed -i \'/# Third party databases/a\\' + ss + '\' ' + pro_path + '/SNPEFF.config\n')
f.write('java -jar ${snpeff} build -gff3 -v ' + '.'.join(GENOME.split('.')[0:-1]) + ' -c ' + pro_path + '/SNPEFF.config\n')
f.write('java -Xmx16g -jar ${snpeff} eff -v ' + '.'.join(GENOME.split('.')[0:-1]) + ' -c ' + pro_path + '/SNPEFF.config -i vcf ' + path_04 + 'all.gatk.vcf >' + path_05 + 'annote.vcf\n')
f.close()
### 以上代码的功能是VCF注释
f = open('./bin/08_run_cnvnator.sh','w')
f.write(script_header1 + '\n')
f.write('for i in')
for i in list3:
    f.write(' ' + i)
f.write('\ndo\n')
f.write('$cnvnator -root ${i}\'.root\' -genome ' + GPATH + GENOME + ' -tree ' + path_03 + '${i}.rmdup.bam\n')
f.write('$cnvnator -root ${i}\'.root\' -genome ' + GPATH + GENOME + ' -his 50 -d ' + REF_DIR + '\n')
f.write('$cnvnator -root ${i}\'.root\' -genome ' + GPATH + GENOME + ' -stat 50\n')
f.write('$cnvnator -root ${i}\'.root\' -genome ' + GPATH + GENOME + ' -partition 50\n')
f.write('$cnvnator -root ${i}\'.root\' -genome ' + GPATH + GENOME + ' -call 50 >' + path_07 + '${i}.cnv\n')
f.write('done\n')
f.close()
#获取运行cnvnator的shell脚本

f = open('./bin/09_depth_cov.sh','w')
f.write(script_header1 + '\n')
for i in list3:
    f.write('$iTools Fqtools stat -InFq ' + RPATH + i + R1 + " -InFq " + RPATH + i + R2 + ' -OutStat ' + RPATH + i + '.baseq\n')
    f.write('$iTools Fqtools stat -InFq ' + CPATH + i + C1 + " -InFq " + CPATH + i + C2 + ' -OutStat ' + CPATH + i + '.baseq\n')
    f.write('$samtools depth -aa ' + path_03 + i + '.rmdup.bam' + ' -o ' + path_03 + i + '.depth\n')
    f.write('awk -F"\\t" \'BEGIN{a=0;b=0;c=0 } {a+=1;b+=$3;if($3 != "0"){c+=1}} END{printf"sample\\tcoverage\\tdepth\\n' + i + '\\t%s\\t%s\\n",c/a,b/a}\' ' + path_03 + i + '.depth >' + path_03 + i + '.depth.stat\n')
f.close()

## stat depth and coverage
###***********************************
f = open('./bin/10_sep_snp_indel.py','w')
pyscript='''#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
import re
pat = re.compile(r'^[ACGT]$')


f = open(sys.argv[1],'r')
f1 = open("all.snp.vcf",'w')
f2 = open('all.indel.vcf','w')



for i in f:
    if i.startswith("##"):
        continue
    if '/2' in i:
        continue
    if i.startswith("#"):
        f1.write(i)
        f2.write(i)
    else:
        temp = i.strip().split()
        if re.search(pat,temp[3]) and re.search(pat,temp[4]):
            f1.write(i)
        else:
            f2.write(i)

f.close()
f1.close()
f2.close()
'''

f.write(pyscript)
f.close()

f = open('./bin/11.get_gene_id.py','w')
py_scripy_1 = r'''#!/data/user/niexiaoqing/anaconda3/bin/python
import sys
import re

def usage():
    print("第一个参数(sys.argv[1]) : .gff文件，读取\n第二个参数(sys.argv[2]) : .cnv文件，读取\n第三个参数(sys.argv[2]) : 指定输出文件的路径，写操作\n注意：该脚本的代码可能会根据.gff文件的数据结构的差异而修改")

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
            gene_list.add(gff_list[a1][3])
            if a1 == gff_len -1:
                break
            a1 += 1

f3 = open(sys.argv[3],'w')
for i in gene_list:
    f3.write(i + "\n")
f3.close()
'''
f.write(py_scripy_1)
f.close()


a='''#!/bin/bash
./01_filter.sh
./02_run_fastqc.sh
./03_get_index.sh   
./04_run_bwa.sh 
./05_run_samtools.sh
./06_run_gatk.sh      
./07_run_snpeff.sh   
./08_run_cnvnator.sh
./09_depth_cov.sh
./10_sep_snp_indel.py {my_gatk}
for i in {my_samples}
do
./11.get_gene_id.py {my_gff} {path_07}$i.cnv $i.cnv.annotation.txt
done

mv snpEff* {my_path_8}
rm *root 
mv *annotation.txt (my_path_8)
'''.format(my_gatk = path_04 + 'all.gatk.vcf',my_path_8 = path_08,path_07 = path_07,my_samples = " ".join(list3),my_gff = gff_path + '/genes.gff')

f = open('./bin/RUN_all.sh','w')
f.write(a)
f.close()

p1 = subprocess.Popen("chmod 755 " + pro_path + '/bin/*',shell=True)
p1.wait()
