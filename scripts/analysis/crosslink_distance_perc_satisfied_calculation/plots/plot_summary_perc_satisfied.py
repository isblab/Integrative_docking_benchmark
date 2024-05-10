import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

def get_perc_across_models(file_path, flag):
    total_perc = []

    with open(file_path, 'r') as file:
        for line in file:
            models_perc = float(line.split()[1])
            total_perc.append(models_perc)

    max_xlink_perc = max(total_perc)
    models_with_max_xlink_sat = sum(1 for perc in total_perc if perc >= max_xlink_perc)
    total_lines = len(total_perc)

    avg_perc = np.mean(total_perc)
    # print(avg_perc, file_path)
    # print(max_xlink_perc if max_xlink_perc < 75 else None)

    if flag == 'max':
        return max_xlink_perc

def get_xlinks_satisfied_and_corres_mdls(name, flag):
    file1_path = os.path.join('/home/muskaan/easal/imp_output/xl_satisfaction/', name + '_perc_satisfied.txt')
    file2_path = os.path.join('/home/muskaan/easal/easal_output/xl_satisfaction/', name + '_perc_satisfied.txt')
    try:
        num = int(name.split('_')[2])
    except:
        num = int(name.split('_')[3])

    # Maximum crosslink satisfaction percentage
    if flag == 'max':
        max_xlink_perc_imp = get_perc_across_models(file1_path, flag) #For imp models
        max_xlink_perc_easal = get_perc_across_models(file2_path, flag) #For easal models
        return max_xlink_perc_imp, max_xlink_perc_easal

    # elif flag == 'avg':
        # avg_perc_imp = get_perc_across_models(file1_path, flag) #For imp models
        # avg_perc_easal = get_perc_across_models(file2_path, flag) #For easal models
        # return avg_perc_imp, avg_perc_easal




input_cases = [["1clv_DSSO_2", "1dfj_DSSO_3", "1r0r_DSSO_3", "1kxp_DSSO_4", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5"],
    ["1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "1dfj_DSSO_9","2b42_DSSO_10", "2hle_DSSO_10"],
    [ "1kxp_DSSO_11", "1dfj_DSSO_12", "2ayo_DSSO_13", "2hle_DSSO_14"],
    ["1dfj_EDC_4", "2ayo_EDC_5", "1r0r_EDC_6", "1kxp_EDC_7", "1clv_EDC_8", "2hle_EDC_9","2b42_EDC_10"],
    ["roca_putc_DSSO_2", "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "gcvpa_gcvpb_DSSO_5", "phes_phet_DSSO_8"],
    ["1dfj_DSSO_12", "2b42_EDC_10","2hle_DSSO_5","phes_phet_DSSO_8"] ]

flag = 'max'
fig, ax = plt.subplots(figsize=(8, 8))
colors = ['#c1d11f','#6ec007', '#00610e', 'red', 'blue']
#
legend_elements = [Line2D([0], [0], color='#c1d11f', label='<5 simulated (D) crosslinks'),
                   Line2D([0], [0], color='#6ec007', label='6-10 simulated (D) crosslinks'),
                   Line2D([0], [0], color='#00610e', label='>10 simulated (D) crosslinks'),
                   Line2D([0], [0], color='red', label='Simulated (E) crosslinks '),
                   Line2D([0], [0], color='blue', label='Experimental (D) crosslinks')]

for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
    for idx, case in enumerate(ic):
        max_xlink_perc_imp, max_xlink_perc_easal = get_xlinks_satisfied_and_corres_mdls(case, flag)
        if max_xlink_perc_imp == 100 and max_xlink_perc_easal == 100:
            plt.scatter(max_xlink_perc_imp, max_xlink_perc_easal, color = color, s =150)
        else:
            plt.scatter(max_xlink_perc_imp, max_xlink_perc_easal, color = color)
        # plt.scatter(max_xlink_perc_easal, color = 'orange')
        print(max_xlink_perc_imp, max_xlink_perc_easal, case)

plt.xlabel('Highest percentage of crosslinks\n satisfied by an IMP model (%)', fontsize=16)
plt.ylabel('Highest percentage of crosslinks\n satisfied by an EASAL model (%)',fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(0, 105)
plt.ylim(0, 105)
plt.legend(handles=legend_elements, fontsize=14)
plt.savefig('/home/muskaan/easal/plots/summary/F2.xlink_per_sat_max.png', dpi=600)
# plt.show()

# flag = 'avg'
# fig, ax = plt.subplots(figsize=(8, 8))
# colors = ['#c1d11f','#6ec007', '#00610e', 'red', 'blue']
#
# legend_elements = [Line2D([0], [0], color='#c1d11f', label='<5 simulated (D) crosslinks'),
#                    Line2D([0], [0], color='#6ec007', label='6-10 simulated (D) crosslinks'),
#                    Line2D([0], [0], color='#00610e', label='>10 simulated (D) crosslinks'),
#                    Line2D([0], [0], color='red', label='Simulated (E) crosslinks '),
#                    Line2D([0], [0], color='blue', label='Experimental (D) crosslinks')]
#
# for color_idx, (ic, color) in enumerate(zip(input_cases, colors)):
#     for idx, case in enumerate(ic):
#         avg_imp, avg_easal = get_xlinks_satisfied_and_corres_mdls(case, flag)
#         # print(avg_imp, avg_easal, case)
#         plt.scatter(avg_imp, avg_easal, color=color)
#
# plt.xlabel('Average percentage of crosslinks\n satisfied in IMP ensemble (%)', fontsize=16)
# plt.ylabel('Average percentage of crosslinks\n satisfied in EASAL ensemble (%)', fontsize=16)
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.xlim(0, 105)
# plt.ylim(0, 105)
# plt.legend(handles=legend_elements, fontsize=14)
# plt.savefig('/home/muskaan/easal/plots/summary/F2.xlink_per_sat_avg.png', dpi=600)
# # plt.show()
