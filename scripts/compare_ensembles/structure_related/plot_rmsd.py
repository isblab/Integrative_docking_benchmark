import matplotlib.pyplot as plt
import os
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

def read_file_and_get_rmsd(name, flag):
    mdl = os.path.join('~/easal/plots/structure_related/rmsd/', name + '.txt')

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
                all_rmsd_imp = [float(x.strip(',[]')) for x in line.split()[2:]]  # Slice from index 3 to end and remove comma
                within_10A_count_imp = sum(1 for rmsd in all_rmsd_imp if abs(rmsd - min_rmsd_imp) <= 5)
                # print(within_10A_count_imp, len(all_rmsd_imp))

            if 'rmsd easal' in line:
                all_rmsd_easal = [float(x.strip(',[]')) for x in line.split()[2:]]  # Slice from index 3 to end
                within_10A_count_easal = sum(1 for rmsd in all_rmsd_easal if abs(rmsd - min_rmsd_easal) <= 5)
                # print(within_10A_count_easal, len(all_rmsd_easal))

    if flag == 'min':
        return min_rmsd_imp, min_rmsd_easal

    if flag == 'all':
        return all_rmsd_imp, all_rmsd_easal

    if flag == 'avg':
        return np.mean(all_rmsd_imp), np.mean(all_rmsd_easal)



#All input cases
input_cases = [["1clv_DSSO_2", "1dfj_DSSO_3", "1r0r_DSSO_3", "1kxp_DSSO_4", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "1dfj_DSSO_9","2b42_DSSO_10", "2hle_DSSO_10"],
    [ "1kxp_DSSO_11", "1dfj_DSSO_12", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_DMTMM_4", "2ayo_DMTMM_5", "1r0r_DMTMM_6", "1kxp_DMTMM_7", "1clv_DMTMM_8", "2hle_DMTMM_9","2b42_DMTMM_10"],
    ["roca_putc_DSSO_2", "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "gcvpa_gcvpb_DSSO_5", "phes_phet_DSSO_8"]]

#Input cases for complexwise plots

#Plotting

colors = ['#c1d11f','#6ec007', '#00610e', 'red', 'purple']

legend_elements = [Line2D([0], [0], color='#c1d11f', label='<5 simulated (DS) crosslinks'),
                   Line2D([0], [0], color='#6ec007', label='6-10 simulated (DS) crosslinks'),
                   Line2D([0], [0], color='#00610e', label='>10 simulated (DS) crosslinks'),
                   Line2D([0], [0], color='red', label='Simulated (DM) crosslinks '),
                   Line2D([0], [0], color='purple', label='Experimental (DS) crosslinks')]

flag = 'min'
count = 0
fig, ax = plt.subplots(figsize=(8, 8))
for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
    for idx, case in enumerate(ic):
        rmsd_imp, rmsd_easal = read_file_and_get_rmsd(case, flag)
        print(case, rmsd_imp, rmsd_easal)

        # if rmsd_imp < 10:
        #     count +=1
        plt.scatter(rmsd_imp, rmsd_easal, color=color)

    # print(count)
    plt.plot([0, 100], [0, 100], '--', color='gray')
    plt.xlabel('Minimum RMSD in IMP ensemble (Å)',fontsize=16)
    plt.ylabel('Minimum RMSD in EASAL ensemble (Å)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlim(0, 80)
    plt.ylim(0, 80)
    plt.legend(handles=legend_elements, fontsize=14)
    plt.savefig('~/easal/plots/structure_related/F5.minimum_rmsd.png',dpi=600)
    # plt.show()

flag = 'all'

for idx, cases in enumerate(input_cases):
    fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})
    for case_idx, case in enumerate(cases):
        row = case_idx // 3
        col = case_idx % 3
        rmsd_imp, rmsd_easal = read_file_and_get_rmsd(case, flag)
        if case.count('_') == 3:  # Checking if the case has 3 underscores
            title = case.split('_')[0] + '_' +case.split('_')[1] + '/' +'/'.join(case.split('_')[2:])
        elif case.count('_') == 2:
            title = '/'.join(case.split('_'))
        axs[row,col].violinplot(rmsd_imp, showmeans=False, showmedians=False)
        axs[row,col].violinplot(rmsd_easal, showmeans=False, showmedians=False)
        axs[row,col].set_title(f'{title}', fontsize=20)
        axs[row,col].set_xlabel('Density',fontsize=18)
        axs[row,col].set_ylabel('RMSD (Å)',fontsize=18)
        axs[row, col].set_ylim(0, 140)
        axs[row,col].tick_params(axis='both', which='major', labelsize=16)
        axs[row, col].legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])
    for i in range(len(cases), 9):
        fig.delaxes(axs.flatten()[i])
    plt.savefig(f'~/easal/plots/structure_related/F5.{idx}.png')
    # plt.show()

flag = 'avg'
fig, ax = plt.subplots(figsize=(10, 8))
double_avg_imp, double_avg_easal = [],[]

for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
    for idx, case in enumerate(ic):
        avg_imp, avg_easal = read_file_and_get_rmsd(case, flag)
        double_avg_imp.append(avg_imp)
        double_avg_easal.append(avg_easal)
        # print(case, avg_imp, avg_easal)

        plt.scatter(avg_imp, avg_easal, color=color)

print(np.mean(double_avg_imp), np.mean(double_avg_easal))
plt.plot([0, 100], [0, 100], '--', color='gray')
plt.xlabel('Average RMSD in IMP ensemble (Å)',fontsize=16)
plt.ylabel('Average RMSD in EASAL ensemble (Å)',fontsize=16)
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(handles=legend_elements, fontsize=14)
plt.savefig('~/easal/plots/structure_related/F5.avg_rmsd.png', dpi=600)
# plt.show()
