# %%
file = '/Users/tomoyauchiyama/code/SD0566/all_alignments_SD0566_01_a_hifiasm_consensus.tsv'

#S1      E1      S2      E2      Reference       Contig  IDY     Ambiguous       Best_group
#2110099 2113376 344613  347896  chr4_1_4500000  ptg000001l      95.16           True
#CONTIG  ptg000002l      567678  mis_unaligned

ctgID = {}
with open(file, 'r') as ctg:
    for line in ctg.readlines():
        line = line.rstrip('\n')
        tmp = line.split()
        if (tmp[0] == 'CONTIG') & (tmp[-1] == 'correct_unaligned'): #correct
            ctgID[tmp[1]] = '-'

outf = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_2.txt'
cnt = {}
with open(outf, 'w') as wf:
    with open(file, 'r') as ctg:
        for line in ctg.readlines():
            line = line.rstrip('\n')
            tmp = line.split()
            if str.isdigit(tmp[0]):
                if ctgID.get(tmp[5]):
                    if cnt.get(tmp[5]):
                        cnt[tmp[5]] += 1
                    else:
                        cnt[tmp[5]] = 1

                    if int(tmp[2]) < int(tmp[3]):
                        strand = '+'
                    else:
                        strand = '-'
                    ctg_name = tmp[5] + '_' + str(cnt[tmp[5]])
                    column = tmp[4] + '\t' + tmp[0] + '\t' + tmp[1] + '\t' + ctg_name + '\t' + tmp[6] + '\t' + strand + '\t' + tmp[2] + '\t' + tmp[3] + '\t' + tmp[7]
                    wf.write(f'{column}')
                    length = int(tmp[1]) - int(tmp[0]) + 1
                    wf.write(f'\t{length}\n')

# %%
"""

sort -n -k 2,2 SD0566_correct_contigs_1.txt > SD0566_correct_contigs_1_sorted.txt
sort -n -k 2,2 SD0566_correct_contigs_2.txt > SD0566_correct_contigs_2_sorted.txt

"""

# %%
file = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_2_sorted.txt'
#chr4_1_4500000  151742  152709  ptg004498l_1    121127  120160  98.76   True    968
#chr4_1_4500000  151742  152709  ptg005330l_1    482     1448    99.07   True    968
#chr4_1_4500000  14806   15814   ptg003956l_1    95.74   -       414327  413321  True    1009
nr = {}
with open(file, 'r') as ctg:
    for line in ctg.readlines():
        line = line.rstrip('\n')
        tmp = line.split()
        id = tmp[1]
        ln = int(tmp[-1])
        idt = float(tmp[4])
        if nr.get(id):
            if ln < int(tmp[-1]):
                ln = int(tmp[-1])
                idt = float(tmp[4])
                nr[id] = line
            elif ln == int(tmp[-1]):
                if idt <= float(tmp[4]):
                    ln = int(tmp[-1])
                    idt = float(tmp[4])
                    nr[id] = line
        else:
            nr[id] = line

outf = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_2_mod.txt'
with open(outf, 'w') as wf:
    for key in sorted(nr.keys(), key=int):
        line = nr[key]
        wf.write(f'{line}\n')

# %%
"""

bedtools merge -i SD0566_correct_contigs_1_sorted.txt -c 4 -o collapse > SD0566_correct_contigs_1_merged.txt
bedtools merge -i SD0566_correct_contigs_2_mod.txt -S + -c 4 -o collapse > SD0566_correct_contigs_2_merged.txt # '+'の方が処理がしやすい、処理の実行時間も早い

"""

# %%
#chr4_1_4500000  12883   545454  ptg003010l_2,ptg003010l_1,ptg003764l_1

file = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_1_merged.txt'
#chr4_1_4500000  12883   545454  ptg003010l_2,ptg003010l_1,ptg003764l_1
#chr4_1_4500000  561860  566647  ptg006636l_1

full_ctg = {}
perfect_start = []
perfect_end = []
with open(file, 'r') as mgd:
    for line in mgd.readlines():
        line = line.rstrip('\n')
        tmp = line.split()

        id = str(tmp[1]) + ',' + str(tmp[2])
        full_ctg[id] = line

        perfect_start.append(int(tmp[1]))
        perfect_end.append(int(tmp[2]))

# %%
gap_regions = []
for i in range(len(perfect_start)-1):
    gap_region = str(perfect_end[i]) + '-' + str(perfect_start[i + 1])
    gap_regions.append(gap_region)

# %%
print(gap_regions)

# %%
import re

file = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_2_merged.txt'

with open(file, 'r') as mgd:
    for line in mgd.readlines():
        line = line.rstrip('\n')
        tmp = line.split()
        if int(tmp[2]) < int(perfect_start[0]):
            id = str(tmp[1]) + ',' + str(tmp[2])
            full_ctg[id] = line
        elif int(perfect_end[-1]) < int(tmp[1]):
            id = str(tmp[1]) + ',' + str(tmp[2])
            full_ctg[id] = line
        else:
            for gap_reg in gap_regions:
                gap_start = re.compile(r'(.*)-(.*)').search(gap_reg).group(1)
                gap_end = re.compile(r'(.*)-(.*)').search(gap_reg).group(2)
                if (int(gap_start) < int(tmp[1])) & (int(tmp[2]) < int(gap_end)):
                    id = str(tmp[1]) + ',' + str(tmp[2])
                    full_ctg[id] = line


# %%
print(len(sorted(full_ctg.keys())))

# %%

"""

cat SD0566_correct_contigs_1_sorted.txt SD0566_correct_contigs_2_mod.txt > SD0566_correct_contigs_12_sorted.txt

"""

# %%
file = '/Users/tomoyauchiyama/code/SD0566/SD0566_correct_contigs_12_sorted.txt'
#chr4_1_4500000  12883   20010   ptg003010l_2    99.06   -       92254   85120   True    7128

chrID = {}
ctgID = {}
ctgdir = {}
with open(file, 'r') as ctg:
    for line in ctg.readlines():
        line = line.rstrip('\n')
        tmp = line.split()
        id = tmp[3]
        chrID[id] = str(tmp[1]) + '-' + str(tmp[2])
        ctgID[id] = str(tmp[6]) + '-' + str(tmp[7])
        ctgdir[id] = tmp[5]


outf = '/Users/tomoyauchiyama/code/SD0566/SD0566_full_correct_contigs.txt'
with open(outf, 'w') as wf:
    wf.write('#chr_id\tchr_start\tchr_end\tctg_id\tctgs_pos(chr)\tctgs_pos(ctg)\tctgs_strand\n')
    for key in sorted(full_ctg.keys(), key=lambda x:int((re.search(r'(\d+),(\d+)', x)).group(1))):
        line = full_ctg[key]
        tmp = line.split()
        chr_loc = [chrID[i] for i in tmp[-1].split(',')]
        chr_pos = ','.join(chr_loc)
        ctg_loc = [ctgID[i] for i in tmp[-1].split(',')]
        ctg_pos = ','.join(ctg_loc)
        ctg_str = [ctgdir[i] for i in tmp[-1].split(',')]
        ctg_dir = ','.join(ctg_str)
        wf.write(f'{line}\t{chr_pos}\t{ctg_pos}\t{ctg_dir}\t')

# %%
s = '12883-20010,19034-68810,22478-545454'
loc = s.split('-')

for i in range(1, len(loc) - 1):
    print(f'{loc[i]}')

# %%
# contig内で切り取る配列領域をcontig間のoverlap領域を除いたものに修正

file = '/Users/tomoyauchiyama/code/SD0566/SD0566_full_correct_contigs.txt'
# chr4_1_4500000  12883   545454
# ptg003010l_2,ptg003010l_1,ptg003764l_1  12883-20010,19034-68810,22478-545454    92254-85120,79816-30132,516982-1        -,-,-
# chr4_1_4500000  2929486 3275331
# ptg006052l_1,ptg004629l_1       2929486-3126723,3106351-3275331 1-197564,1-168955       +,+

outf = '/Users/tomoyauchiyama/code/SD0566/SD0566_full_correct_contigs_mod.txt'
with open(outf, 'w') as wf:
    wf.write('#chr_id\tchr_start\tchr_end\tctg_id\tctgs_pos(chr)\tctgs_pos(ctg)\tctgs_strand\tmod_tgs_pos(ctg)\n')
    with open(file, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            if re.compile(r'^#').search(line):
                pass
            else:
                tmp = line.split()
                chr_loc = tmp[4].split('-')
                ovrlp_bps = []
                mod_ctg_loc = []

                if len(chr_loc) == 2:
                    if '+' in tmp[-1]:
                        mod_ctg_pos = tmp[5]
                    else:
                        ed = re.compile(r'(\d+)-(\d+)').search(tmp[5]).group(1)
                        st = re.compile(r'(\d+)-(\d+)').search(tmp[5]).group(2)
                        mod_ctg_pos = str(st) + '-' + str(ed)

                elif len(chr_loc) > 2:
                    for i in range(1, len(chr_loc) - 1):
                        ed_i0 = re.compile(r'(\d+),(\d+)').search(chr_loc[i]).group(1)
                        st_i1 = re.compile(r'(\d+),(\d+)').search(chr_loc[i]).group(2)
                        if int(st_i1) <= int(ed_i0):
                            ovrlp_len = int(ed_i0) - int(st_i1) + 1
                            ovrlp_bps.append(ovrlp_len)
                        else:
                            ovrlp_len = 0
                            ovrlp_bps.append(ovrlp_len)

                    ctg_loc = tmp[5].split(',')
                    for i in range(len(ctg_loc)):
                        if i == 0:
                            if '+' in tmp[-1]:
                                mod_loc = ctg_loc[i]
                            else:
                                ed = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(1)
                                st = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(2)
                                mod_loc = str(st) + '-' + str(ed)
                        else:
                            if '+' in tmp[-1]:
                                st = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(1)
                                ed = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(2)
                                mod_st = int(st) + int(ovrlp_bps[i - 1])
                                mod_loc = str(mod_st) + '-' + str(ed)
                            else:
                                #92254-85120,79816-30132,516982-1
                                ed = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(1)
                                st = re.compile(r'(\d+)-(\d+)').search(ctg_loc[i]).group(2)
                                mod_ed = int(ed) - int(ovrlp_bps[i - 1])
                                mod_loc = str(st) + '-' + str(mod_ed)
                        mod_ctg_loc.append(mod_loc)

                    mod_ctg_pos = ','.join(mod_ctg_loc)

                wf.write(f'{line}\t{mod_ctg_pos}\n')

# %%
file = '/Users/tomoyauchiyama/code/SD0566/all_alignments_SD0566_01_a_hifiasm_consensus.tsv'
outf = '/Users/tomoyauchiyama/code/SD0566/SD0566_01_a_hifiasm_consensus_pos.tsv'
with open(outf, 'w') as wf:
    with open(file, 'r') as ctg:
        for line in ctg.readlines():
            line = line.rstrip('\n')
            tmp = line.split()
            if str.isdigit(tmp[0]):
                wf.write(f'{line}\n')
# %%
