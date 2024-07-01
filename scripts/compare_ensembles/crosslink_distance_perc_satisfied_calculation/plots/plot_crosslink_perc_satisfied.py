import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
import argparse

def calc_crosslink_satisfaction(file_path):
    percentages = []
    with open(file_path, 'r') as file:
        for line in file:
            percentage = float(line.split(' ')[1].strip())
            percentages.append(percentage)

    return percentages

def reading_file_and_get_perc(name):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/xl_satisfaction/', name + '_perc_satisfied.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/xl_satisfaction/', name + '_perc_satisfied.txt')

    perc_imp = calc_crosslink_satisfaction(file1_path)
    perc_easal = calc_crosslink_satisfaction(file2_path)

    return perc_imp, perc_easal
#TODO instead of commenting out, use sys.argv to choose what to plot each time. Reduces human errors.
#Input cases
parser = argparse.ArgumentParser(description='Input cases.')

parser.add_argument('--lt5', action='store_true', help='DSSO less than 5')
parser.add_argument('--b6_10', action='store_true', help='DSSO 6-10')
parser.add_argument('--mt10', action='store_true', help='DSSO more than 10')
parser.add_argument('--exp', action='store_true', help='DSSO experimental')
parser.add_argument('--edc', action='store_true', help='EDC')
parser.add_argument('--sel', action='store_true', help='Selected cases')

args = parser.parse_args()
if args.lt5:
    input_cases = ["1clv_DSSO_2", "1dfj_DSSO_3", "1r0r_DSSO_3", "1kxp_DSSO_4", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"]
    outf = 'S1.less_than_5'
elif args.b6_10:
    input_cases = ["1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "1dfj_DSSO_9","2b42_DSSO_10", "2hle_DSSO_10"]
    outf = 'S1.between_6_to_10'
elif args.mt10:
    input_cases = [ "1kxp_DSSO_11", "1dfj_DSSO_12", "2ayo_DSSO_13", "2hle_DSSO_14"]
    outf = 'S1.more_than_10'
elif args.exp:
    input_cases = ["roca_putc_DSSO_2", "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "gcvpa_gcvpb_DSSO_5", "phes_phet_DSSO_8"]
    outf = 'S1.experimental'
elif args.edc:
    input_cases = ["1dfj_EDC_4", "2ayo_EDC_5", "1r0r_EDC_6", "1kxp_EDC_7", "1clv_EDC_8", "2hle_EDC_9","2b42_EDC_10"]
    outf = 'S1.EDC'
elif args.sel:
    input_cases = ["2b42_DSSO_5", "1dfj_DSSO_9", "roca_putc_DSSO_2", "2hle_EDC_9"]
    outf = 'F2.xlink_per_sat_complexwise'

#Plotting

fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})

for idx, case in enumerate(input_cases):
    row = idx // 3
    col = idx % 3

    perc_imp, perc_easal = reading_file_and_get_perc(case)
    if case.count('_') == 3:  # Checking if the case has 3 underscores
        title = case.split('_')[0] + '_' +case.split('_')[1] + '/' +'/'.join(case.split('_')[2:])
    elif case.count('_') == 2:
        title = '/'.join(case.split('_'))
    axs[row, col].violinplot(perc_imp, showmeans=False, showmedians=False, widths = 1)
    axs[row, col].violinplot(perc_easal, showmeans=False, showmedians=False)
    axs[row, col].set_title(f'{title}', fontsize=20)
    axs[row, col].set_ylabel('Percentage of\n crosslinks satisfied (%)', fontsize=18)
    axs[row, col].set_xlabel('Density', fontsize=18)
    axs[row, col].set_ylim(0, 110)
    axs[row, col].tick_params(axis='both', which='major', labelsize=14)
    axs[row, col].legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])

# Remove empty subplots
for i in range(len(input_cases), 9):
    fig.delaxes(axs.flatten()[i])

# plt.savefig(f'/home/muskaan/easal/plots/percentage_satisfied/{outf}.png')
plt.show()
