import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def read_file_and_get_rmsd(name, flag):
    mdl = os.path.join('/home/muskaan/easal/plots/structure_related/rmsd/', name + '.txt')

    with open(mdl, 'r') as f:
        for line in f:
            if flag == 'min':
                if 'min imp' in line:
                    min_rmsd_imp = line.split(' ')[-1]
                elif 'min easal' in line:
                    min_rmsd_easal = line.split(' ')[-1]
                    return float(min_rmsd_imp), float(min_rmsd_easal)

            if flag == 'all':
                if 'rmsd imp' in line:
                    all_rmsd_imp = [float(x.strip(',[]')) for x in line.split()[3:]]  # Slice from index 3 to end and remove comma
                elif 'rmsd easal' in line:
                    all_rmsd_easal = [float(x.strip(',[]')) for x in line.split()[3:]]  # Slice from index 3 to end
                    return all_rmsd_imp, all_rmsd_easal

#All input cases
input_cases = [["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"],
    ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"],
    ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]]

#Input cases for complexwise plots
# selected_cases = ["roca_putc_DSSO_2", "gcvpa_gcvpb_DSSO_5", "1clv_EDC_8", "1dfj_EDC_4", "1clv_DSSO_2", "2b42_DSSO_5", "1dfj_DSSO_9", "2hle_DSSO_14", "1dfj_DSSO_12"]
flag = sys.argv[1] #Specify whether to plot 'summary' or 'complexwise'
#Plotting

if flag == 'min':

    fig, ax = plt.subplots(figsize=(8, 8))

    colors = ['r', 'g', 'b', 'm', 'orange']

    legend_elements = [Line2D([0], [0], color='red', label='DSSO simulated less than 5'),
                       Line2D([0], [0], color='green', label='DSSO simulated 6-10'),
                       Line2D([0], [0], color='blue', label='DSSO simulated more than 10'),
                       Line2D([0], [0], color='magenta', label='EDC simulated'),
                       Line2D([0], [0], color='orange', label='DSSO experimental')]

    for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
        for idx, case in enumerate(ic):
            rmsd_imp, rmsd_easal = read_file_and_get_rmsd(case, flag)
            print(rmsd_imp, rmsd_easal)
            plt.scatter(rmsd_imp, rmsd_easal, color=color)

    plt.xlabel('Minimum RMSD in IMP model (Å)',fontsize=14)
    plt.ylabel('Minimum RMSD in EASAL model (Å)',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlim(0, 60)
    plt.ylim(0, 60)
    plt.legend(handles=legend_elements)
    plt.savefig('/home/muskaan/easal/plots/structure_related/F5.minimum_rmsd.png',dpi=600)
    plt.show()


elif flag == 'all':
    for idx, cases in enumerate(input_cases):
        fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})
        for case_idx, case in enumerate(cases):
            row = case_idx // 3
            col = case_idx % 3
            rmsd_imp, rmsd_easal = read_file_and_get_rmsd(case, flag)

            axs[row,col].violinplot(rmsd_imp, showmeans=False, showmedians=False)
            axs[row,col].violinplot(rmsd_easal, showmeans=False, showmedians=False)
            axs[row,col].set_title(f'{case}')
            axs[row,col].set_xlabel('Density',fontsize=14)
            axs[row,col].set_ylabel('Crosslink distance diff\n in native structure\n vs model (Å)',fontsize=14)
            axs[row,col].tick_params(axis='both', which='major', labelsize=12)
            axs[row, col].legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])
        for i in range(len(cases), 9):
            fig.delaxes(axs.flatten()[i])
        plt.savefig(f'/home/muskaan/easal/plots/structure_related/F5.{idx}.png')
        plt.show()
