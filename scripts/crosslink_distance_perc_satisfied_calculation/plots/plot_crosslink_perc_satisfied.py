import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot_crosslink_satisfaction(ax, file_path, label, color, name):
    percentages = []
    with open(file_path, 'r') as file:
        for line in file:
            percentage = float(line.split(' ')[1])
            percentages.append(percentage)

    ax.hist(percentages, bins=10, color=color, label=label, alpha=0.7, histtype='step',linewidth=2)

    ax.set_title(f'{name}')
    ax.set_xlabel('Percentage')
    ax.set_ylabel('Number of models')
    ax.legend()

def file_parsing(name, ax):
    file1_path = os.path.join('/home/muskaan/easal_imp/xl_satisfaction/', name + '_percentage_satisfied.txt')
    file2_path = os.path.join('/home/muskaan/easal_output_muskaan/xl_satisfaction/', name + '_perc_satisfied.txt')

    plot_crosslink_satisfaction(ax, file1_path, 'IMP', color='blue', name=name)
    plot_crosslink_satisfaction(ax, file2_path, 'EASAL', color='orange', name=name)

# input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"] #DSSO less than 5
# input_cases = ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"] #DSSO 6-10
# input_cases = ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"] #DSSO more than 10
input_cases = ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4"] #DSSO experimental
# input_cases = ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"] #EDC
# all_input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5","1dfj_DSSO_9",
#                     "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
#                     "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
#                      "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4",
#                      "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"]
fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})

for idx, case in enumerate(input_cases):
    row = idx // 3
    col = idx % 3
    file_parsing(case, axs[row, col])

plt.savefig('/home/muskaan/easal_plots/percentage_satisfied/all_subplots.png')
plt.show()
