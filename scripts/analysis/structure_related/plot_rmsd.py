import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def read_file_and_get_rmsd(name, flag):
    mdl = os.path.join('/home/muskaan/easal/plots/structure_related/rmsd/', name + '.txt')

     # Initialize variables to store RMSD values
    min_rmsd_imp = None
    min_rmsd_easal = None
    all_rmsd_imp = None
    all_rmsd_easal = None


    with open(mdl, 'r') as f:
        for line in f:
            if 'min imp' in line:
                min_rmsd_imp = float(line.split(' ')[-1])
            if 'min easal' in line:
                min_rmsd_easal = float(line.split(' ')[-1])

            if 'rmsd imp' in line:
                all_rmsd_imp = [float(x.strip(',[]')) for x in line.split()[3:]]  # Slice from index 3 to end and remove comma
                within_10A_count_imp = sum(1 for rmsd in all_rmsd_imp if abs(rmsd - min_rmsd_imp) <= 10)
                # print(within_10A_count_imp, len(all_rmsd_imp))

            if 'rmsd easal' in line:
                all_rmsd_easal = [float(x.strip(',[]')) for x in line.split()[3:]]  # Slice from index 3 to end
                within_10A_count_easal = sum(1 for rmsd in all_rmsd_easal if abs(rmsd - min_rmsd_easal) <= 10)
                # print(within_10A_count_easal, len(all_rmsd_easal))

    if flag == 'min':
        return min_rmsd_imp, min_rmsd_easal

    if flag == 'all':
        return all_rmsd_imp, all_rmsd_easal

    if flag == 'within_10A':
        return min_rmsd_imp, (within_10A_count_imp/len(all_rmsd_imp))*100, min_rmsd_easal, (within_10A_count_easal/len(all_rmsd_easal))*100



#All input cases
input_cases = [["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"],
    ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"],
    ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"],
    ["1dfj_DSSO_9", "2ayo_EDC_5","1dfj_DSSO_12","phes_phet_DSSO_8"] ]

#Input cases for complexwise plots
flag = sys.argv[1] #Specify whether to plot 'minimum_rmsd' or 'all rmsd'
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
            print(case, rmsd_imp, rmsd_easal)
            plt.scatter(rmsd_imp, rmsd_easal, color=color)

    plt.xlabel('Minimum RMSD in IMP model (Å)',fontsize=14)
    plt.ylabel('Minimum RMSD in EASAL model (Å)',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlim(0, 80)
    plt.ylim(0, 80)
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
            title = '/'.join(case.split('_'))
            axs[row,col].violinplot(rmsd_imp, showmeans=False, showmedians=False)
            axs[row,col].violinplot(rmsd_easal, showmeans=False, showmedians=False)
            axs[row,col].set_title(f'{title}', fontsize=20)
            axs[row,col].set_xlabel('Density',fontsize=18)
            axs[row,col].set_ylabel('RMSD(Å)',fontsize=18)
            axs[row, col].set_ylim(0, 130)
            axs[row,col].tick_params(axis='both', which='major', labelsize=16)
            axs[row, col].legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])
        for i in range(len(cases), 9):
            fig.delaxes(axs.flatten()[i])
        plt.savefig(f'/home/muskaan/easal/plots/structure_related/F5.{idx}.png')
        # plt.show()

elif flag == 'within_10A':

    fig, ax = plt.subplots(figsize=(10, 8))

    for idx, cases in enumerate(input_cases):
        for case_idx, case in enumerate(cases):
            min_rmsd_imp, num_mdls_imp, min_rmsd_easal, num_mdls_easal = read_file_and_get_rmsd(case, flag)
            print(case, min_rmsd_imp, num_mdls_imp, min_rmsd_easal, num_mdls_easal)
            plt.scatter(min_rmsd_imp, num_mdls_imp, color = 'blue', label='IMP')
            plt.scatter(min_rmsd_easal, num_mdls_easal, color = 'orange', label='EASAL')

    plt.xlabel('Minimum RMSD in models (Å)',fontsize=16)
    plt.ylabel('Percentage of models\n within 10Å of minimum_rmsd (%)', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])
    plt.savefig('/home/muskaan/easal/plots/structure_related/F5.models_within_10_min_rmsd.png', dpi=600)
    # plt.show()
