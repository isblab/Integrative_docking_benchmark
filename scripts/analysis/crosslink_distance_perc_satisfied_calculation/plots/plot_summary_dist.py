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

fig, ax = plt.subplots(figsize=(8, 8))
sat_imp, sat_easal = 0,0

for case in input_cases:
    avg_imp, avg_easal = read_file_and_get_dist(case)
    if 'DSSO' in case:
        marker = 'o'
        color = 'red'
        if avg_imp < 32:
            sat_imp+=1
        if avg_easal < 32:
            sat_easal += 1
    elif 'EDC' in case:
        marker = '^'
        color = 'green'
        if avg_imp < 20:
            sat_imp+=1
        if avg_easal < 20:
            sat_easal += 1

    plt.scatter(avg_imp, avg_easal, marker=marker, color=color)
    # print(case, avg_imp, avg_easal)
print(sat_imp, sat_easal) # Number of cases with avg distance within the cutoff

plt.xlabel('Average Crosslink Distance in IMP Ensemble (Å)', fontsize=14)
plt.ylabel('Average Crosslink Distance in EASAL Ensemble (Å)', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(10, 60)
plt.ylim(10, 60)


legend_elements = [Line2D([0], [0], marker='o', color='red', label='DSSO', markersize=10, linestyle='None'),
                   Line2D([0], [0], marker='^', color='green', label='EDC', markersize=10, linestyle='None')]

plt.legend(handles=legend_elements)
plt.savefig('/home/muskaan/easal/plots/F3.xlink_dist_summary.png',dpi=600)
plt.show()
