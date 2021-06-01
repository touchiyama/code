import os
import re
import string as str

#cwd=os.getcwd()
#print(cwd)

dna={}
with open('../Arabidopsis/TAIR10/test.fa','r') as fa:
    for line in fa.readlines():
        line=line.rstrip('\n')
        if re.compile(r'^>').search(line): #染色体番号を取得
            tmp=line.split()
            #m=re.compile(r'^>([0-9]+)').search(tmp[0])
            m=re.compile(r'^>(.*)').search(tmp[0])
            chr='Chr'+ m.group(1)
        else: #ゲノム配列を取得する
            if dna.get(chr):
                dna[chr]+=line
                #print(chr)
            else:
                dna[chr]=line

"""
for k, v in dna.items(): #keyとvalueを同時に表示させる
    print(k, v)
"""

#Chr1    TAIR10  CDS     3760    3913    .       +       0       Parent=AT1G01010.1,AT1G01010.1-Protein;
CDS={}
with open('../Arabidopsis/TAIR10/TAIR10_GFF3_genes.gff','r') as gff:
    for line in gff.readlines():
        line=line.rstrip('\n')
        tmp=line.split()
        if tmp[2] == 'CDS':
            ele=tmp[8].split(',')
            m=re.compile(r'Parent=(.*)').search(ele[0])
            id=m.group(1)+','+tmp[6]
            loc=tmp[0]+':'+tmp[3]+'-'+tmp[4]
            if CDS.get(id):
                CDS[id]=CDS[id]+'#'+loc
            else:
                CDS[id]=loc

cdn={'TTT' : 'F' , 'TTC' : 'F' , 'TTA' : 'L' , 'TTG' : 'L' , 'CTT' : 'L' , 'CTC' : 'L' , 'CTA' : 'L' , 'CTG' : 'L' ,
        'ATT' : 'I' , 'ATC' : 'I' , 'ATA' : 'I' , 'ATG' : 'M' , 'GTT' : 'V' , 'GTC' : 'V' , 'GTA' : 'V' , 'GTG' : 'V' ,
        'TCT' : 'S' , 'TCC' : 'S' , 'TCA' : 'S' , 'TCG' : 'S' , 'CCT' : 'P' , 'CCC' : 'P' , 'CCA' : 'P' , 'CCG' : 'P' ,
        'ACT' : 'T' , 'ACC' : 'T', 'ACA' : 'T' , 'ACG' : 'T' , 'GCT' : 'A' , 'GCC' : 'A' , 'GCA' : 'A' , 'GCG' : 'A' ,
        'TAT' : 'Y' , 'TAC' : 'Y', 'TAA' : '*' , 'TAG' : '*' , 'CAT' : 'H' , 'CAC' : 'H' , 'CAA' : 'Q' , 'CAG' : 'Q' ,
        'AAT' : 'N' , 'AAC' : 'N' , 'AAA' : 'K' , 'AAG' : 'K' , 'GAT' : 'D' , 'GAC' : 'D' , 'GAA' : 'E' , 'GAG' : 'E',
        'TGT' : 'C' , 'TGC' : 'C' , 'TGA' : '*' , 'TGG' : 'W' , 'CGT' : 'R' , 'CGC' : 'R' , 'CGA' : 'R' , 'CGG' : 'R' ,
        'AGT' : 'S' , 'AGC' : 'S' , 'AGA' : 'R' , 'AGG' : 'R' , 'GGT' : 'G' , 'GGC' : 'G' , 'GGA' : 'G' , 'GGG' : 'G' }

for k, v in CDS.items():
    print('>'+k)
    ele=v.split('#')
    seq=[]
    for i in ele:
        m=re.compile(r'(.*):(\d+)-(\d+)').search(i)
        chr=m.group(1)
        start=int(m.group(2))-1
        end=int(m.group(3))
        dna=dna[chr]
        seq.append(dna[start:end])

    tmp=k.split(',')
    if tmp[1] == '+':
        nucl=''.join(seq)
    else:
        nucl=''.join(seq.reverse())
        nucl=reversed(seq)
        nucl=seq.translate(str.maketrans('ATGC', 'TACG'))

    for i in range(0,int(len(nucl)/3)):
        triple=nucl[3*i:3*i+3]
        if(triple in cdn):
            print(cdn[triple], end="")
        else:
            print ("Z")
    print('\n')






