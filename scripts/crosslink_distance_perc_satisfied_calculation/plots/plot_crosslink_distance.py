import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
import argparse

def calc_crosslink_distances(ax, file_path, label, color, name):
    distances = []
    with open(file_path, 'r') as file:
        for line in file:
            distance = float(line.split()[1])
            distances.append(distance)

    return distances

#TODO see comments for the plot_crosslink_perc_satisfied and follow the same

def read_file_and_get_dist(name, ax):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/crosslink_distances/', name + '_distances.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', name + '_distances.txt')

    dist_imp = calc_crosslink_distances(ax, file1_path, 'IMP', color='blue', name=name)
    dist_easal = calc_crosslink_distances(ax, file2_path, 'EASAL', color='orange',name=name)

    return dist_imp, dist_easal

#input cases
parser = argparse.ArgumentParser(description='Input cases.')

parser.add_argument('--less_than_5', help='DSSO less than 5')
parser.add_argument('--between_6_to_10', help='DSSO 6-10')
parser.add_argument('--more_than_10', help='DSSO more than 10')
parser.add_argument('--experimental', help='DSSO experimental')
parser.add_argument('--edc', help='EDC')
parser.add_argument('--selected', help='Selected cases')

args = parser.parse_args()
if args.less_than_5:
    input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"]
elif args.between_6_to_10:
    input_cases = ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"]
elif args.more_than_10:
    input_cases = ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"]
elif args.experimental:
    input_cases = ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5", "roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]
elif args.edc:
    input_cases = ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"]
elif args.selected:
    input_cases = ["1dfj_DSSO_9", "1clv_DSSO_2", "1dfj_EDC_4", "gcvpa_gcvpb_DSSO_5"]

#Plotting
fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})

for idx, case in enumerate(input_cases):
    row = idx // 3
    col = idx % 3
    dist_imp, dist_easal = read_file_and_get_dist(case, axs[row, col])

    axs[row,col].violinplot(dist_imp, showmeans=False, showmedians=False)
    axs[row,col].violinplot(dist_easal, showmeans=False, showmedians=False)
    axs[row,col].set_title(f'{case}')
    axs[row,col].set_xlabel('Density',fontsize=14)
    axs[row,col].set_ylabel('Distance between\n crosslinked residues\n in model (Å)',fontsize=14)
    axs[row,col].tick_params(axis='both', which='major', labelsize=12)
    axs[row, col].legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])

plt.savefig(f'/home/muskaan/easal/plots/distance_distribution/{sys.argv[1]}.png', dpi=600)
plt.show()
