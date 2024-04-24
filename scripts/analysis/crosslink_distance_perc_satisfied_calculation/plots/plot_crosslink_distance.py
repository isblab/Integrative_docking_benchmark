import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
import argparse

def calc_crosslink_distances(file_path):
    distances = []
    with open(file_path, 'r') as file:
        for line in file:
            distance = float(line.split()[1].strip())
            distances.append(distance)

    return distances

#TODO see comments for the plot_crosslink_perc_satisfied and follow the same

def read_file_and_get_dist(name):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/crosslink_distances/', name + '_distances.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', name + '_distances.txt')

    dist_imp = calc_crosslink_distances(file1_path)
    dist_easal = calc_crosslink_distances(file2_path)

    return dist_imp, dist_easal

#input cases
parser = argparse.ArgumentParser(description='Input cases.')

parser.add_argument('--lt5', action='store_true', help='DSSO less than 5')
parser.add_argument('--b6_10', action='store_true', help='DSSO 6-10')
parser.add_argument('--mt10', action='store_true', help='DSSO more than 10')
parser.add_argument('--exp', action='store_true', help='DSSO experimental')
parser.add_argument('--edc', action='store_true', help='EDC')
parser.add_argument('--sel', action='store_true', help='Selected cases')

args = parser.parse_args()
if args.lt5:
    input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"]
    outf = 'S2.less_than_5'
elif args.b6_10:
    input_cases = ["1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10"]
    outf = 'S2.between_6_to_10'
elif args.mt10:
    input_cases = ["1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14"]
    outf = 'S2.more_than_10'
elif args.exp:
    input_cases = ["gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5", "roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]
    outf = 'S2.experimental'
elif args.edc:
    input_cases = ["1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9"]
    outf = 'S2.EDC'
elif args.sel:
    input_cases = ["gcvpa_gcvpb_DSSO_5","1clv_DSSO_2","1dfj_EDC_4", "1dfj_DSSO_9"]
    outf = 'F3.xlink_dist_complexwise'

#Plotting

fig, axs = plt.subplots(2, 2, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})

for idx, case in enumerate(input_cases):
    row = idx // 2
    col = idx % 2

    dist_imp, dist_easal = read_file_and_get_dist(case)
    title = '/'.join(case.split('_'))
    axs[row, col].violinplot(dist_imp, showmeans=False, showmedians=False)
    axs[row, col].violinplot(dist_easal, showmeans=False, showmedians=False)
    axs[row, col].set_title(f'{title}', fontsize=24)
    axs[row, col].set_xlabel('Density (A.U.)', fontsize=20)
    axs[row, col].set_ylabel('Crosslink\n distances(Å)', fontsize=20)
    axs[row, col].set_ylim(0, 120)
    axs[row, col].tick_params(axis='both', which='major', labelsize=16)
    axs[row, col].legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])

# Remove empty subplots
#for i in range(len(input_cases), 9):
 #   fig.delaxes(axs.flatten()[i])


plt.savefig(f'/home/muskaan/easal/plots/distance_distribution/violin/{outf}.png')
# plt.show()
