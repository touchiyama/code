file = '/Users/tomoyauchiyama/code/all_alignments_SD0566_01_a_hifiasm_consensus.tsv'

#S1      E1      S2      E2      Reference       Contig  IDY     Ambiguous       Best_group
#2110099 2113376 344613  347896  chr4_1_4500000  ptg000001l      95.16           True
#CONTIG  ptg000002l      567678  mis_unaligned
ctgID = {}
with open(file, 'r') as ctg:
    for line in ctg.readlines():
        line = line.rstrip('\n')
        tmp = line.split()
        if (tmp[0] == 'CONTIG') & (tmp[-1] == 'correct'):
            ctgID[tmp[1]] = '-'

outf = './SD0566_correct_contigs_2.txt'
with open(outf, 'w') as wf:
    with open(file, 'r') as ctg:
        for line in ctg.readlines():
            line = line.rstrip('\n')
            tmp = line.split()
            if str.isdigit(tmp[0]):
                if ctgID.get(tmp[5]):
                    wf.write(f'{line}')
                    length = int(tmp[1]) - int(tmp[0]) + 1
                    wf.write(f'\t{length}\n')

sort -n -k 1,1 SD0566_correct_contigs_2.txt > SD0566_correct_contigs_2_sorted.txt
SD0566_correct_contigs_2_sorted.txtの内容を処理して、１つのsccafoldを作る。
具体的に以下の処理を行う。
(1) 間隙領域には、Nの文字を入れる
(2) alignment領域に、それに対応するcontig配列から切り出してきた塩基を挿入
(3) (2) でend > startの表示になっている時、逆相補鎖に変換
(4) 完璧に近づけようと思うならば、間隙領域の補間を行う。
