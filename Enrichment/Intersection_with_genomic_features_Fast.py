coord_basename = '/EPIC_hg19/EPIC/Processed_files/'
coord_dirs = ['Dnase_HS_and_open_chromatin', 'Histones_HepG2_(Broad)', 'Histones_HepG2_(UW)', 'Others','TFBS_HepG2_(HAIB_and_SYDH)']
input_file = open(r'/../Meth_map_goruped.txt', 'r')
output_dir = '/../Intersections_with_genomic_features/'

import os
from operator import itemgetter
import copy

def mergeIntervals(all_intervals):
    sorted_intervals = sorted(all_intervals, key = itemgetter(0, 1, 2))
    merged_intervals = []
    for interval in sorted_intervals:
        if not merged_intervals:
            merged_intervals.append(interval)
        else:
            previous = merged_intervals[-1]
            prev_chr = previous[0]
            prev_start = previous[1]
            prev_end = previous[2]
            cur_chr = interval[0]
            cur_start = interval[1]
            cur_end = interval[2]
            if prev_chr == cur_chr and cur_start <= (prev_end):   # overlap detected!
                # (prev_end + 1) -> non overlapping, but adjacent intervals are also considered as overlapping
                # (prev_end) -> only intervals overlapping for at least 1 bp are considered;
                min_start = min(prev_start, cur_start)
                max_end = max(prev_end, cur_end)
                merged_intervals[-1] = [cur_chr, min_start, max_end]
            else:
                merged_intervals.append(interval)
    return merged_intervals


##reading input file
header, data = [], []
for line in input_file:
    if line[-1] == '\n':
        line = line[:-1]
    L = line.split('\t')
    if line[0] == '#':
        header = [L[0], L[1], L[2]]# ['#chr', 'coord', 'CpGSite']  its header
    else:
        data.append([L[0], int(L[1]), L[2]])# [['chr1', 15865, 'cg13869341'], ...]  its info, all lines that are not header
input_file.close()
if not header:
    header = ['#col1', 'col2', 'col3']
print('Input file loaded;')

##sorting data 
data = sorted(data, key=itemgetter(0, 1))
print('Input data sorted;')


for coord_dir in coord_dirs:
    curr_header = copy.copy(header)
    curr_data = copy.deepcopy(data)
    print('\n')
    print(coord_dir)
    print('curr_header:', curr_header)#############
    print('curr_data[1]:', curr_data[1])#############

    ##output file
    output_file = open(output_dir + 'target_vs ' + coord_dir + '.txt', 'w')
    workdir = coord_basename + coord_dir +"/"
    filenames = [f for f in os.listdir(workdir) if not f.startswith('.')]
    for filename in filenames:# for each genomic feature
        curr_header.append(filename[:-4])
        coord_file = open(workdir + filename, 'r')
        print(coord_file)
        intervals = []
        stats = [0, 0, 0, 0]# Yes, 1Kb, 2Kb, No
        for line in coord_file:
            if line[0] == '#':
                continue
            if line[-1] == '\n':
                line = line[:-1]
            L = line.split('\t')
            print(L)
            intervals.append([L[0], int(L[1]), int(L[2])])

        coord_file.close()
        print('\t', filename[:-4], '\t', len(intervals), 'intervals loaded;')
        merged_intervals = mergeIntervals(intervals)
        if len(merged_intervals) < len(intervals):
            print('\t', len(merged_intervals), 'non-overlapping intervals;')
        i, j = 0, 0
        while i < len(curr_data):
            curr_cpg = curr_data[i]
            if j < len(merged_intervals):
                curr_interval = merged_intervals[j]
                if (curr_cpg[0] < curr_interval[0]) or (curr_cpg[0] == curr_interval[0] and curr_cpg[1] < curr_interval[1]-2000):# CpG is located before the interval start
                    curr_data[i].append('No')
                    stats[3] += 1
                    i += 1# skip current cpg and continue with the next cpg
                elif curr_cpg[0] == curr_interval[0] and curr_cpg[1] >= curr_interval[1]-2000 and curr_cpg[1] <= curr_interval[2]+2000:# hit!
                    if curr_cpg[1] >= curr_interval[1] and curr_cpg[1] <= curr_interval[2]:
                        curr_data[i].append('Yes')
                        stats[0] += 1
                    else:# current CpG is within 2Kb zone, but there might be a better match with the next interval...
                        k = 0
                        all_decisions = []
                        while (j+k) < len(merged_intervals):
                            next_interval = merged_intervals[j+k]
                            decision = False
                            if curr_cpg[0] == next_interval[0]:
                                if curr_cpg[1] >= next_interval[1] and curr_cpg[1] <= next_interval[2]:
                                    decision = 0
                                elif curr_cpg[1] >= next_interval[1]-1000 and curr_cpg[1] <= next_interval[2]+1000:
                                    decision = 1
                                elif curr_cpg[1] >= next_interval[1]-2000 and curr_cpg[1] <= next_interval[2]+2000:
                                    decision = 2
                            if decision:
                                all_decisions.append(decision)
                                k += 1
                            else:
                                break
                        final_decision = min(all_decisions)
                        if final_decision == 1:
                            curr_data[i].append('1Kb')
                            stats[1] += 1
                        elif final_decision == 2:
                            curr_data[i].append('2Kb')
                            stats[2] += 1
                        elif final_decision == 0:
                            curr_data[i].append('Yes')
                            stats[0] += 1
                    i += 1# continue with the next cpg
                else:# CpG is located either after the interval end, or on the next chromosome
                    j += 1# skip interval and continue with the next one
            else:
                break
        # process the remaining CpGs, if any
        while i < len(curr_data):
            curr_data[i].append('No')
            stats[3] += 1
            i += 1
        print('\t', stats)
        """
        troubleshoot = []
        value, count = None, 0
        for cpg in curr_data:
            length = len(cpg)
            if value:
                if length == value:
                    count += 1
                else:
                    troubleshoot.append([value, count])
                    value = length
                    count = 1
            else:
                value = length
                count += 1
        troubleshoot.append([value, count])
        print('\t', troubleshoot)
        if len(troubleshoot) > 1:
            for item in troubleshoot[1:]:
                short = item[0]
                for i in range(len(curr_data)):
                    cpg = curr_data[i]
                    if len(cpg) == short:
                        print('\t\t', i, cpg[:3])
        """
    header_line = ''
    for field in curr_header:
        header_line = header_line + field + '\t'
    output_file.write(header_line[:-1] + '\n')
    for cpg in curr_data:
        output_line = ''
        for field in cpg:
            output_line = output_line + str(field) + '\t'
        output_file.write(output_line[:-1] + '\n')
    output_file.flush()
    output_file.close()
    del curr_header
    del curr_data
    print('\tOutput file written;')
print('Done!')
