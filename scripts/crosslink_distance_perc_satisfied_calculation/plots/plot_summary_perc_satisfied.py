import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def get_perc_across_models(file_path, num):
    perc_more_than_80 = 0

    with open(file_path, 'r') as file:
        for line in file:
            models_perc = float(line.split()[1])
            if models_perc > 80.0:
                perc_more_than_80 += 1
    total_lines = sum(1 for line in open(file_path))
    # print((perc_more_than_80/total_lines) * 100)
    return (perc_more_than_80/total_lines) * 100

def file_parsing(name):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/xl_satisfaction/', name + '_perc_satisfied.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/xl_satisfaction/', name + '_perc_satisfied.txt')

    try:
        num = int(name.split('_')[2])
        labels = name.split('_')[0]+name.split('_')[2]
    except:
        num = int(name.split('_')[3])
        labels = name.split('_')[0]+'_'+name.split('_')[1]+name.split('_')[3]

    avg_imp = get_perc_across_models(file1_path, num)
    avg_easal = get_perc_across_models(file2_path, num)

    if 'DSSO' in name:
        marker = 'o'
        color = 'red'
    elif 'EDC' in name:
        marker = '^'
        color = 'green'

    plt.scatter(avg_imp, avg_easal, marker=marker, color=color)
    # plt.annotate(labels, (avg_imp, avg_easal), textcoords="offset points", xytext=(10,-10), ha='center')

input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4"]

for case in input_cases:
    file_parsing(case)

# plt.rcParams['font.family'] = 'Arial'
plt.xlabel('% models in IMP Ensemble satisfying more than 80% crosslinks')
plt.ylabel('% models in EASAL Ensemble satisfying more than 80% crosslinks')

legend_elements = [Line2D([0], [0], marker='o', color='red', label='DSSO', markersize=10, linestyle='None'),
                   Line2D([0], [0], marker='^', color='green', label='EDC', markersize=10, linestyle='None')]

plt.legend(handles=legend_elements)
plt.savefig('/home/muskaan/easal/plots/summary_plot_perc.png')
plt.show()
