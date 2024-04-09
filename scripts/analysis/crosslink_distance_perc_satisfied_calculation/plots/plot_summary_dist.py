import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def get_avg_dist(file_path, num):
    distances = []

    with open(file_path, 'r') as file:
        for line in file:
            distance = float(line.split()[1])
            distances.append(distance)

    avg_distances = []
    for i in range(0, len(distances), num):
        avg_distance = np.mean(distances[i:i+num])
        avg_distances.append(avg_distance) #This is avg distance across all crosslinks

    return np.mean(avg_distances) #This is avg distance across all models

def read_file_and_get_dist(name):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/crosslink_distances/', name + '_distances.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', name + '_distances.txt')

    try:
        num = int(name.split('_')[2])
    except:
        num = int(name.split('_')[3])

    avg_imp = get_avg_dist(file1_path, num)
    avg_easal = get_avg_dist(file2_path, num)

    return avg_imp, avg_easal


input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]

for case in input_cases:
    #TODO make the function name more descriptive of what it is doing than file parsing.
    #TODO again do all plotting in 1 function so that the plt object is in one place.

    avg_imp, avg_easal = read_file_and_get_dist(case)

    if 'DSSO' in case:
        marker = 'o'
        color = 'red'
    elif 'EDC' in case:
        marker = '^'
        color = 'green'

    plt.scatter(avg_imp, avg_easal, marker=marker, color=color)

# plt.rcParams['font.family'] = 'Arial'
plt.xlabel('Average Crosslink Distance in IMP Ensemble (Å)')
plt.ylabel('Average Crosslink Distance in EASAL Ensemble (Å)')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(10, 50)
plt.ylim(10, 50)


legend_elements = [Line2D([0], [0], marker='o', color='red', label='DSSO', markersize=10, linestyle='None'),
                   Line2D([0], [0], marker='^', color='green', label='EDC', markersize=10, linestyle='None')]

plt.legend(handles=legend_elements)
plt.savefig('/home/muskaan/easal/plots/summary_plot_dist.png',dpi=600)
plt.show()
