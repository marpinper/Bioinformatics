feature_dir = '/../Intersections_with_genomic_features/'
feature_cpg_columns = [1, 2, 3]# chr1, 245524538, cg03130248
decisions = ['Yes']# positive flags ('Yes', '1Kb', '2Kb')
input_file = open(r'/../hmC_up.txt', 'r')
input_cpg_columns = [5, 6, 1]# chr1, 245524538, cg03130248   
filtering = 'float(L[10-1]) > 0'# expression or empty string
output_file = open(r'/../hmC_up_vs_genomic_features.txt', 'w')

import os
from scipy.stats import fisher_exact

cpg_subset = []
line_count = 0
for line in input_file:
    if 'CpG' in line:
        continue
    if line[-1] == '\n':
        line = line[:-1]
    L = line.split('\t')
    print(L)
    cpg_id = ''
    line_count += 1
    for column in input_cpg_columns:
        print(column)
        print(L[column-1])
        cpg_id = cpg_id + L[column-1] + '|'
    cpg_id = cpg_id[:-1]
    print(cpg_id)

    print(filtering)
    print(eval(filtering))
    print(L[10])
    print(L[1])
    print(L[10-1])
    if filtering and not eval(filtering):

        continue
    cpg_subset.append(cpg_id)
input_file.close()
print(line_count, 'input lines processed;')
print(len(cpg_subset), 'input CpGs loaded;')

output_file.write('#UCSC_feature\tHit_450K\tN_450K\t%450K\tHit_subset\tN_subset\t%subset\tEnrichment_R\tFisher_pval\n')
filenames = [f for f in os.listdir(feature_dir) if not f.startswith('.')]
for filename in filenames:
    print(filename[:-4])
    output_file.write('\n')
    output_file.write('Feature file = ' + filename + ':\n')
    feature_file = open(feature_dir + filename, 'r')
    all_data, relevant_data, header = {}, {}, ''
    for line in feature_file:
        if line[-1] == '\n':
            line = line[:-1]
        L = line.split('\t')
        if line[0] == '#':
            header = L[3:]
            continue
        cpg_id = ''
        for column in feature_cpg_columns:
            cpg_id = cpg_id + L[column-1] + '|'
        cpg_id = cpg_id[:-1]
        all_data[cpg_id] = L[3:]
    feature_file.close()
    print('Data loaded;')
    not_found = 0
    for cpg in cpg_subset:
        if cpg in all_data:
            relevant_data[cpg] = all_data[cpg]
        else:
            not_found += 1
    print(len(relevant_data), 'relevant data extracted;')
    if not_found > 0:
        print(not_found, 'relevant CpGs not found!')
    stats = []
    if not header:
        print('No header!')
        continue
    for i in range(len(header)):
        curr_feature = header[i]
        print('\t' + curr_feature)
        pos_all, pos_rel = 0, 0
        for key in all_data:
            flag = all_data[key][i]
            if flag in decisions[0]:# this is a positive flag
                pos_all += 1
        for key in relevant_data:
            flag = relevant_data[key][i]
            if flag in decisions[0]:
                pos_rel += 1
        if pos_all == 0 or pos_rel == 0:
            enrichment = 'NA'
        else:
            enrichment = round((pos_rel/len(relevant_data)) / (pos_all/len(all_data)), 3)
        oddsratio, pvalue = fisher_exact([[pos_rel, len(relevant_data)-pos_rel], [pos_all, len(all_data)-pos_all]])
        stats.append([curr_feature, pos_all, len(all_data), round(pos_all/len(all_data)*100, 2), pos_rel, len(relevant_data), round(pos_rel/(len(relevant_data)+0.001)*100, 2), enrichment, pvalue])
    #stats = sorted(stats, key = lambda item: abs(1-item[7]), reverse=True)
    for record in stats:
        output_line = ''
        for field in record:
            output_line = output_line + str(field) + '\t'
        output_file.write(output_line[:-1] + '\n')
output_file.flush()
output_file.close()
print('Done!')
