import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def get_avg_dist(file_path, num, fp):
    distances = []

    with open(file_path, 'r') as file:
        for line in file:
            distance = float(line.split()[1])
            distances.append(distance)

    per_model = []
    for i in range(0, len(distances), num):
        fp_dist = distances[i:i+num][-fp:] # FP are in the end of the xlinks csv file 
        sat = 0
        for j in range(len(fp_dist)):
            if fp_dist[j] < 32.0: # Check if distance is less than violation threshold or not
                sat += 1
        per_model.append((sat/fp)*100)

    return np.mean(per_model)

def read_file_and_get_dist(name):
    file1_path = os.path.join('~/easal/imp_output/test_fp/crosslink_distances/', name + '_distances.txt')
    file2_path = os.path.join('~/easal/easal_output/test_fp/crosslink_distances/', name + '_distances.txt')

    num = int(name.split('_')[2])
    fp = 1 if num == 5 else 2 if num == 10 else 3

    avg_imp = get_avg_dist(file1_path, num, fp)
    avg_easal = get_avg_dist(file2_path, num, fp)

    return avg_imp, avg_easal


input_cases = ['1dfj_DSSO_10', '1clv_DSSO_5', '1dfj_DSSO_11', '1dfj_DSSO_10', '1kxp_DSSO_10', '1r0r_DSSO_5', '2ayo_DSSO_15', '2b42_DSSO_10', '2hle_DSSO_10', '2hle_DSSO_15', '2ayo_DSSO_10']

fig, ax = plt.subplots(figsize=(8, 8))
sat_imp, sat_easal = 0,0

for case in input_cases:
    avg_imp, avg_easal = read_file_and_get_dist(case)
    plt.scatter(avg_imp, avg_easal, color = 'green')
    print(case, avg_imp, avg_easal)

plt.plot([0,110], [0, 110], '--', color='gray')
plt.xlabel('Average FP Crosslink satisfication in IMP Ensemble (Å)', fontsize=16)
plt.ylabel('Average FP Crosslink satisfication in wall-EASAL Ensemble (Å)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(0, 110)
plt.ylim(0, 110)
plt.savefig('~/easal/plots/summary/S6.fp_xlink.png',dpi=600)
plt.show()
