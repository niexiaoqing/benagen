#!/data/user/niexiaoqing/anaconda3/bin/python
import re
import pandas as pd
import sys
#name = sys.argv[1]
f = open(sys.argv[1],'r')
str1 = f.read()
f.close()

pat = re.compile(r'Number variants by type.*?(<table.*?</table>)',re.S)
pat1 = re.compile(r'Number of effects by impact.*?(<table.*?</table>)',re.S)
pat2 = re.compile(r'Number of effects by functional class.*?(<table.*?</table>)',re.S)

table1_1 = re.search(pat,str1).group(1)
table1_2 = pd.read_html(table1_1)[0]
table1_3 = table1_2.drop(index=(table1_2.loc[(table1_2['Total']==0)].index))
table1_3["Total"] = table1_3["Total"].map(lambda x:format(x,","))
table1_3.to_csv("snp_stat.txt",index=False,sep="\t")
print(table1_3)
print("\n\n")

table2_1 = re.search(pat1,str1).group(1)
table2_2 = pd.read_html(table2_1)[0]
table2_2["Count"] = table2_2["Count"].map(lambda x:format(x,","))
table2_2 = table2_2.drop(['Unnamed: 1'], axis=1)
table2_2.rename(columns={'Type (alphabetical order)':'Type'},inplace=True)
table2_2.to_csv("affect_var.txt",sep="\t",index=False)
print(table2_2)
print("\n\n")


table3_1 = re.search(pat2,str1).group(1)
table3_2 = pd.read_html(table3_1)[0]
table3_2["Count"] = table3_2["Count"].map(lambda x:format(x,","))
table3_2 = table3_2.drop(['Unnamed: 1'], axis=1)
table3_2.rename(columns={'Type (alphabetical order)':'Type'},inplace=True)
table3_2.to_csv("var_class.txt",index=False,sep="\t")
print(table3_2)
print("\n\n")
