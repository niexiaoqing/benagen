#!/data/user/niexiaoqing/anaconda3/bin/python
import subprocess
import re
import pandas as pd
import sys

def usage():
    if len(sys.argv) != 2:
        print("""USAGE:
get_report.py xxxx.clean.fastqc.html
### the only args specify a html file from which we get a image""") 
        exit(1)

usage()



pat = re.compile(r'#ReadNum:[^0-9]*(\d*)[^0-9]*BaseNum')
pat1 = re.compile(r'BaseNum:[^0-9]*(\d*)[^0-9]*ReadLeng:')
pat2 = re.compile(r'ReadLeng:[^0-9]*(\d*)')
pat3 = re.compile(r'#BaseQ:10--20 :.*?>Q20:[^0-9]*(.*)')

p1 = subprocess.Popen('[[ -e outcome ]] || mkdir outcome',shell=True)
p1.wait()
p2 = subprocess.Popen('ls ../01.rawdata/*baseq',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
p2.wait()
str1 = str(p2.stdout.read(),encoding="utf8")
list1 = str1.split()
samples = []

p3 = subprocess.Popen('cp /data/user/niexiaoqing/bena/addition_files/info_model ./outcome/info',shell=True,stdout=subprocess.PIPE)
p3.wait()

p4 = subprocess.Popen('cp /data/user/niexiaoqing/bena/addition_files/hand_write_for_one_project.txt ./',shell=True,stdout=subprocess.PIPE)
p4.wait()

for i in list1:
    temp = i.split('/')[-1]
    temp1 = '.'.join(temp.split('.')[0:-1])
    samples.append(temp1)
    print(temp1)

samp_dic = {}

def get_seq_stat(): 
    global samp_dic
    f = open('./outcome/seq.stat','w')
    f.write("Sample\tRaw paired reads\tRaw bases\tClean paired reads\tClean bases\tQ20\n")
    samp_dic = {}

    for i in samples:
        temp_r = open("../01.rawdata/" + i + '.baseq','r')
        temp_c = open('../02.cleandata/' + i +'.baseq','r')
        temp_r1 = temp_r.read()
        temp_c1 = temp_c.read()
        temp_r.close()
        temp_c.close()
        rdn_r = re.search(pat,temp_r1)
        bsn_r = re.search(pat1,temp_r1)
        rd_len_r = re.search(pat2,temp_r1)
        rdn_c = re.search(pat,temp_c1)
        bsn_c = re.search(pat1,temp_c1)
        rd_len_c = re.search(pat2,temp_c1)
        c_q20 = re.findall(pat3,temp_c1)
        ave = 0.5*(float(c_q20[0].strip('%')) + float(c_q20[1].strip('%')))
        ave = round(ave,2)
        f.write(i + "\t")
        f.write(format(int(rdn_r.group(1)),',') + "\t")
        f.write(format(int(bsn_r.group(1))*2,',') + "\t")
        f.write(format(int(rdn_c.group(1)),',') + "\t")
        f.write(format(int(bsn_c.group(1))*2,',') + "\t")
        f.write(str(ave) + '%\n')
        samp_dic[i] = rdn_r.group(1)
    f.close()    

def get_map_result():
    f = open('./outcome/map.result','w')
    f.write("Sample name\tClean paired reads\tOverall alignment rate\tDepth(x)\tCoverage rate\n")
    pat = re.compile(r'\((.*)% : N/A\)')
    for i in samples:
        temp = open('../log/' + i + '.map_ratio.log','r')
        temp_m = temp.read()
        temp.close()
        temp = open('../03.mapping/' + i + '.depth.stat','r')
        temp_d = temp.readlines()
        depth = temp_d[1].strip().split("\t")[2]
        coverage = temp_d[1].strip().split("\t")[1]
        temp.close()
        print(depth)
        print(coverage)
        map_rate = re.search(pat,temp_m)
        map_r = map_rate.group(1) + '%'
        print(map_r)
        f.write(i + '\t' + format(int(samp_dic[i]),',') + '\t' + map_r + '\t' + '%.02f'%float(depth) + '\t' + '%.03f'%(float(coverage)*100) + '%\n' )
    f.close()


def get_variation_result():
    f = open('./outcome/variation.result','w')
    f.write("Chrom\tPos\tRef\tAlt")
    lines = 0
    f1 = open("../04.snp_indel/all.gatk.vcf",'r')
    pat = re.compile(r'./.:(\d*),(\d*):(\d*)')
    for i in f1:
        if i.startswith("##"):
            continue
        else:
            if i.startswith("#"):
                temp_a = i.strip().split()
                temp_a1 = temp_a[9:]
                for j in temp_a1:
                    f.write("\t" + j)
                f.write("\n")
            else:
                temp_b = i.strip().split()
                temp_b1 = temp_b[9:]
                if '/2' not in i and './.' not in i:
                    if lines == 5:
                        break
                    lines += 1
                    f.write((temp_b[0]) + '\t')
                    f.write((temp_b[1]) + '\t')
                    f.write((temp_b[3]) + '\t')
                    f.write((temp_b[4]))
                    print("*"*80)
                    for jj in temp_b1:
                        temp_b2 = re.search(pat,jj)
                        ref = temp_b2.group(1)
                        alt = temp_b2.group(2)
                        str1 = ref + '_' + alt
                        f.write("\t" + str1)
                    f.write("\n")
    f.close()


def get_all_filter_snpeff():
    varient_type = ["missense_variant","synonymous_variant","synonymous_variant","synonymous_variant","frameshift_variant"]
    f = open('./outcome/all.filter.snpeff','w')
    f.write("Chrom\tPos\tRef\tAlt\tMutation_position_type\tGeneid\tProtein\n")
    f1 = open("../08.annotation/snp_indel/all.snpeff.txt")
    for i in f1:
        if i.startswith("#"):
            continue 
        temp = i.strip().split()
        temp_1 = temp[7].split('|')
        if temp_1[1] in varient_type:
            varient_type.remove(temp_1[1])
            print(temp_1[1])
            print(temp_1[3])
            print(temp_1[10])
            tt = [temp[0],temp[1],temp[3],temp[4],temp_1[1],temp_1[4],temp_1[10]]
            f.write("\t".join(tt) + "\n")
            if varient_type == []:
                break
    f.close()


def get_cnv_result():
    f = open('./outcome/cnv.result','w')
    f.write("Variation_type\tPos\tCNV_length\tE-value\n")
    f1 = open("../07.cnv/" + samples[0] + '.cnv','r')
    lines = 0
    for i in f1:
        temp = i.strip().split()
        temp1 = [temp[0],temp[1],temp[2],temp[4]]
        f.write("\t".join(temp1) + '\n')
        lines +=1
        if lines == 5:
            break
    f.close()


def get_tables():
    f = open("../08.annotation/snp_indel/snpEff_summary.html",'r')
    str1 = f.read()
    f.close()
    
    pat = re.compile(r'Number variants by type.*?(<table.*?</table>)',re.S)
    pat1 = re.compile(r'Number of effects by impact.*?(<table.*?</table>)',re.S)
    pat2 = re.compile(r'Number of effects by functional class.*?(<table.*?</table>)',re.S)
    
    table1_1 = re.search(pat,str1).group(1)
    table1_2 = pd.read_html(table1_1)[0]
    table1_3 = table1_2.drop(index=(table1_2.loc[(table1_2['Total']==0)].index))
    table1_3["Total"] = table1_3["Total"].map(lambda x:format(x,","))
    table1_3.to_csv("./outcome/snp_stat.txt",index=False,sep="\t")
    print(table1_3)
    print("\n\n")
    table2_1 = re.search(pat1,str1).group(1)
    table2_2 = pd.read_html(table2_1)[0]
    table2_2["Count"] = table2_2["Count"].map(lambda x:format(x,","))
    table2_2 = table2_2.drop(['Unnamed: 1'], axis=1)
    table2_2.rename(columns={'Type (alphabetical order)':'Type'},inplace=True)
    table2_2.to_csv("./outcome/affect_var.txt",sep="\t",index=False)
    print(table2_2)
    print("\n\n")
    table3_1 = re.search(pat2,str1).group(1)
    table3_2 = pd.read_html(table3_1)[0]
    table3_2["Count"] = table3_2["Count"].map(lambda x:format(x,","))
    table3_2 = table3_2.drop(['Unnamed: 1'], axis=1)
    table3_2.rename(columns={'Type (alphabetical order)':'Type'},inplace=True)
    table3_2.to_csv("./outcome/var_class.txt",index=False,sep="\t")
    print(table3_2)
    print("\n\n")
        
str2 = ''
def get_img_and_reads_num():
    global str2
    f =open(sys.argv[1],'r',encoding="utf8")
    str1 = f.read()
    f.close()
    pat = re.compile(r'<img class="indented" src="data:image/png;base64,.*?/>',re.S)
    mmy = re.search(pat,str1)
    str2 = mmy.group(0)
    str2 = str2.replace('class="indented"','alt="300x200"')
    str2 = str2.replace('alt="Per base quality graph" width="800" height="600"','')
    print(str2)
    f = open('./outcome/seq.stat','r')
    list1 = f.readlines()[1:]
    f.close()
    r_reads = [int(i.split("\t")[2].replace(",",'')) for i in list1]
    c_reads = [int(i.split("\t")[4].replace(",",'')) for i in list1]
    n1 = sum(r_reads)
    n2 = sum(c_reads)
    n1_1 = round(n1/1000000000,2)
    n2_1 = round(n2/1000000000,2)
    print(n1_1)
    print(n2_1)
    f1 = open('./outcome/project','w')
    f1.write("my_img: \n")
    f1.write(" " + str2 + '\n')
    f1.write("Rreads: \n")
    f1.write(" " + str(n1_1) + "\n")
    f1.write("Creads: \n")
    f1.write(" " + str(n2_1) + "\n")
    f1.write("reference: \n")
    f1.write(" 这是基因组\n")
    f1.write("samples: \n")
    f1.write(" 这是samples")
    f1.close()





get_seq_stat()
get_map_result()
get_variation_result()
get_all_filter_snpeff()
get_cnv_result()
get_tables()
get_img_and_reads_num()
    
    
  
    
