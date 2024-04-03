import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def get_perc_across_models(file_path):
    max_perc = []

    with open(file_path, 'r') as file:
        for line in file:
            models_perc = float(line.split()[1])
            max_perc.append(models_perc)

    max_xlink_perc = max(max_perc)
    models = sum(1 for perc in max_perc if perc >= max_xlink_perc)
    total_lines = len(max_perc)
    return max_xlink_perc, (models / total_lines) * 100


def get_max_xlinks_satisfied_and_corres_mdls(name):
    #TODO more descriptive function name than file parsing: try something likw: get_max_xlinks_satisfied_and_corresponding_models
    #TODO weird to call the plot functions for the same figure in 2 functions, because plt is usually an object.
    #TODO just return numbers from this and do all plotting in main function

    file1_path = os.path.join('/home/muskaan/easal/imp_output/xl_satisfaction/', name + '_perc_satisfied.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/xl_satisfaction/', name + '_perc_satisfied.txt')

    try:
        num = int(name.split('_')[2])
    except:
        num = int(name.split('_')[3])

    #For imp models
    max_xlink_perc_imp, max_models_imp = get_perc_across_models(file1_path)

    #For easal models
    max_xlink_perc_easal, max_models_easal = get_perc_across_models(file2_path)

    return max_xlink_perc_imp, max_models_imp, max_xlink_perc_easal, max_models_easal

fig, ax = plt.subplots(figsize=(8, 8))

input_cases = ["1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]

for case in input_cases:
    max_xlink_perc_imp, max_models_imp, max_xlink_perc_easal, max_models_easal = get_max_xlinks_satisfied_and_corres_mdls(case)
    plt.scatter(max_xlink_perc_imp, max_models_imp, color = 'blue')
    plt.scatter(max_xlink_perc_easal, max_models_easal, color = 'orange')

plt.xlabel('Percentage of maximum crosslinks\nsatisfied by any model', fontsize=12)
plt.ylabel('Percentage of models at\nmaximum crosslink percentage',fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, 105)
plt.ylim(0, 105)
plt.legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])
plt.savefig('/home/muskaan/easal/plots/summary_plot_perc.png', dpi=600)
plt.show()
