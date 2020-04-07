from __future__ import division
# This script is useful for long linear plots, which look not so eye-candy in Excel


input_file = open(r'/../Joined_results_BS_mC.txt', 'r')
data_columns = [ (2, 'blue'), (3, 'black'),(4, 'red'), (5, 'green')]#,(34, 'magenta'), (35, 'pink'), (36, 'brown'), (37, 'gray'),(38, 'yellow'), (39, 'cyan'),(40, 'orchid')]
name_column = 1# int or False
chr_column = False# int or False
output_filename = '/../Joined_results_regions_BS_mC_up_down.png'
figure_title = 'Enrichment rates of hmC CpGs (Wilcoxon) in UCSC genomic features'

#matplotlib inline
import matplotlib.pyplot as plt

for column in data_columns:
    vars()['data_' + str(column[0])] = []
names = []

for line in input_file:
    if line[-1] == '\n':
        line = line[:-1]
    L = line.split('\t')
    if line[0] == '#':
        for column in data_columns:
        
            print(L)

            
            vars()['header_' + str(column[0])] = L[column[0] - 1]
            #vars()['header_' + str(column[0])] = L[column[0] - 1].split('_')[0]
        continue
    for column in data_columns:
        try:
            float(L[column[0] - 1])
        except:
            value = None
        else:
            value = float(L[column[0] - 1])
        vars()['data_' + str(column[0])].append(value)
    if name_column and not chr_column:
        name = L[name_column - 1]
    elif name_column and chr_column:
        name = L[name_column - 1] + ' (' + L[chr_column - 1] + ')'
    names.append(name)
input_file.close()

for column in data_columns:
    current_data = vars()['data_' + str(column[0])]
    plt.plot(range(len(current_data)), current_data, '-o', linewidth = 2, color = column[1], label = vars()['header_' + str(column[0])])
if name_column:
    plt.xticks(range(len(names)), names, rotation = 'vertical', fontsize = 15)
fig = plt.gcf()
fig.set_size_inches(len(current_data)/3, 10)
plt.legend()
plt.title(figure_title)
plt.xlim(xmax = 153, xmin = 0.0)
plt.plot([0, len(current_data)], [1, 1], 'k-')###
plt.savefig(output_filename, bbox_inches='tight', pad_inches=1)
plt.clf()
print('Done!')
