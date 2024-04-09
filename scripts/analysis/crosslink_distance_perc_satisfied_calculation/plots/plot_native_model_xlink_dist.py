import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.lines import Line2D

def get_xlink_dist_diff(mdl, native, num, flag):
    mdl_dist = []
    nat_dist = []
    with open(mdl, 'r') as mdl, open(native, 'r') as nat:
        for line1 in mdl:
            dist1 = float(line1.split()[1])
            mdl_dist.append(dist1)

        for line2 in nat:
            dist2 = float(line2.split()[1])
            nat_dist.append(dist2)

    diff_dist = []

    for i in range(len(mdl_dist)):
        nat_index = i % len(nat_dist)  # Calculate the index in nat_dist
        diff = mdl_dist[i] - nat_dist[nat_index]
        diff_dist.append(abs(diff))

    if flag == 'complexwise':
        return diff_dist
    else:
        return np.mean(diff_dist)

def read_file_and_get_dist(name, flag):
    imp_mdl = os.path.join('/home/muskaan/easal/imp_output/crosslink_distances/', name + '_distances.txt')
    easal_mdl = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', name + '_distances.txt')
    native = os.path.join('/home/muskaan/easal/native_pdb_distances/', name + '_true_structure_distances.txt')

    try:
        num = int(name.split('_')[2]) # Number of crosslinks
    except:
        num = int(name.split('_')[3]) # Number of crosslinks in experimental cases

    dist_easal = get_xlink_dist_diff(easal_mdl, native, num, flag)
    dist_imp = get_xlink_dist_diff(imp_mdl, native, num, flag)

    return dist_imp, dist_easal

#All input cases
input_cases = [["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"],
    ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"],
    ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]]

#Input cases for complexwise plots
selected_cases = ["roca_putc_DSSO_2", "gcvpa_gcvpb_DSSO_5", "1clv_EDC_8", "1dfj_EDC_4", "1clv_DSSO_2", "2b42_DSSO_5", "1dfj_DSSO_9", "2hle_DSSO_14", "1dfj_DSSO_12"]
flag = sys.argv[1] #Specify whether to plot 'summary' or 'complexwise'
#Plotting

if flag == 'summary':

    fig, ax = plt.subplots(figsize=(8, 8))

    colors = ['r', 'g', 'b', 'm', 'orange']

    legend_elements = [Line2D([0], [0], color='red', label='DSSO simulated less than 5'),
                       Line2D([0], [0], color='green', label='DSSO simulated 6-10'),
                       Line2D([0], [0], color='blue', label='DSSO simulated more than 10'),
                       Line2D([0], [0], color='magenta', label='EDC simulated'),
                       Line2D([0], [0], color='orange', label='DSSO experimental')]

    for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
        for idx, case in enumerate(ic):
            dist_imp, dist_easal = read_file_and_get_dist(case, flag)
            plt.scatter(dist_imp, dist_easal, color=color)

    plt.xlabel('Crosslink distance difference in\n the native structure vs IMP model (Å)',fontsize=14)
    plt.ylabel('Crosslink distance difference in\n the native structure vs EASAL model (Å)',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlim(0, 50)
    plt.ylim(0, 50)
    plt.legend(handles=legend_elements)
    plt.savefig('/home/muskaan/easal/plots/summary/F4.xlink_dsit_to_native_summary.png',dpi=600)
    plt.show()

elif flag == 'complexwise':

    fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})
    for idx, case in enumerate(selected_cases):
        row = idx // 3
        col = idx % 3
        dist_imp, dist_easal = read_file_and_get_dist(case, flag)

        axs[row,col].violinplot(dist_imp, showmeans=False, showmedians=False)
        axs[row,col].violinplot(dist_easal, showmeans=False, showmedians=False)
        axs[row,col].set_title(f'{case}')
        axs[row,col].set_xlabel('Density',fontsize=14)
        axs[row,col].set_ylabel('Crosslink distance diff\n in native structure\n vs model (Å)',fontsize=14)
        axs[row,col].tick_params(axis='both', which='major', labelsize=12)
        axs[row, col].legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])
    plt.savefig('/home/muskaan/easal/plots/compare_xlink_dist_native/F4.xlink_dsit_to_native_complexwise.png')
    plt.show()
