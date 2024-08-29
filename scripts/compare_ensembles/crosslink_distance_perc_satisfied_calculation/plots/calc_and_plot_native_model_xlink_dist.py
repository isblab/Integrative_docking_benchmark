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
input_cases = [["1clv_DSSO_2", "1dfj_DSSO_3", "1r0r_DSSO_3", "1kxp_DSSO_4", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "1dfj_DSSO_9","2b42_DSSO_10", "2hle_DSSO_10"],
    [ "1kxp_DSSO_11", "1dfj_DSSO_12", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_DMTMM_4", "2ayo_DMTMM_5", "1r0r_DMTMM_6", "1kxp_DMTMM_7", "1clv_DMTMM_8", "2hle_DMTMM_9","2b42_DMTMM_10"],
    ["roca_putc_DSSO_2", "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "gcvpa_gcvpb_DSSO_5", "phes_phet_DSSO_8"],
    ["1clv_DMTMM_8","2b42_DSSO_5","2hle_DSSO_14", "roca_putc_DSSO_2"]]

#Input cases for complexwise plots
# selected_cases = [["2b42_DSSO_5", "roca_putc_DSSO_2" , "2hle_DSSO_14", "1dfj_DMTMM_4"]]
flag = sys.argv[1] #Specify whether to plot 'summary' or 'complexwise'
#Plotting

if flag == 'summary':

    fig, ax = plt.subplots(figsize=(8, 8))

    colors = ['#c1d11f','#6ec007', '#00610e', 'red', 'purple']

    legend_elements = [Line2D([0], [0], color='#c1d11f', label='<5 simulated (D) crosslinks'),
                       Line2D([0], [0], color='#6ec007', label='6-10 simulated (D) crosslinks'),
                       Line2D([0], [0], color='#00610e', label='>10 simulated (D) crosslinks'),
                       Line2D([0], [0], color='red', label='Simulated (E) crosslinks '),
                       Line2D([0], [0], color='purple', label='Experimental (D) crosslinks')]

    for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
        for idx, case in enumerate(ic):
            dist_imp, dist_easal = read_file_and_get_dist(case, flag)
            plt.scatter(dist_imp, dist_easal, color=color)
            print(case, dist_imp, dist_easal)
    plt.plot([0, 50], [0, 50], '--', color='gray')
    plt.xlabel('Average crosslink distance difference, IMP vs native (Å)',fontsize=16)
    plt.ylabel('Average crosslink distance difference, EASAL vs native (Å)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlim(0, 50)
    plt.ylim(0, 50)
    plt.legend(handles=legend_elements, fontsize = 14)
    plt.savefig('/home/muskaan/easal/plots/summary/F4.xlink_dist_to_native_summary.png',dpi=600)
    # plt.show()

elif flag == 'complexwise':
    for idx, cases in enumerate(input_cases):
        fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})
        for case_idx, case in enumerate(cases):
            row = case_idx // 3
            col = case_idx % 3
            dist_imp, dist_easal = read_file_and_get_dist(case, flag)
            if case.count('_') == 3:  # Checking if the case has 3 underscores
                title = case.split('_')[0] + '_' +case.split('_')[1] + '/' +'/'.join(case.split('_')[2:])
            elif case.count('_') == 2:
                title = '/'.join(case.split('_'))
            axs[row,col].violinplot(dist_imp, showmeans=False, showmedians=False)
            axs[row,col].violinplot(dist_easal, showmeans=False, showmedians=False)
            axs[row,col].set_title(f'{title}', fontsize=20)
            axs[row,col].set_xlabel('Density',fontsize=18)
            axs[row, col].set_ylim(0, 105)
            axs[row,col].set_ylabel('Crosslink distance difference (Å)',fontsize=18)
            axs[row,col].tick_params(axis='both', which='major', labelsize=16)
            axs[row, col].legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])
        for i in range(len(cases), 9):
            fig.delaxes(axs.flatten()[i])
        plt.savefig(f'/home/muskaan/easal/plots/compare_xlink_dist_native/S3.xlink_dist_to_native_{idx}.png')
