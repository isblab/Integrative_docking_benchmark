import matplotlib.pyplot as plt
import os
import sys
import argparse
import matplotlib.patches as mpatches

def get_xlink_dist_diff(mdl, native, num):
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
        # print(nat_dist[nat_index],  mdl_dist[i], diff_dist)
        # exit()

    return diff_dist

#TODO see comments for the plot_crosslink_perc_satisfied and follow the same

def read_file_and_get_dist(name):
    imp_mdl = os.path.join('/home/muskaan/easal/imp_output/crosslink_distances/', name + '_distances.txt')
    easal_mdl = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', name + '_distances.txt')
    native = os.path.join('/home/muskaan/easal/native_pdb_distances/', name + '_true_structure_distances.txt')

    try:
        num = int(name.split('_')[2]) # Number of crosslinks
    except:
        num = int(name.split('_')[3]) # Number of crosslinks in experimental cases

    dist_easal = get_xlink_dist_diff(easal_mdl, native, num)
    dist_imp = get_xlink_dist_diff(imp_mdl, native, num)


    return dist_imp, dist_easal

#input cases
input_cases = ['1clv_DSSO_2', '1dfj_DSSO_9', '1dfj_EDC_4', '2hle_EDC_9', '2b42_DSSO_5', 'roca_putc_DSSO_2', 'gcvpa_gcvpb_DSSO_5']

#Plotting
fig, axs = plt.subplots(3, 3, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.5})

for idx, case in enumerate(input_cases):
    row = idx // 3
    col = idx % 3
    dist_imp, dist_easal = read_file_and_get_dist(case)

    axs[row,col].violinplot(dist_imp, showmeans=False, showmedians=False)
    axs[row,col].violinplot(dist_easal, showmeans=False, showmedians=False)
    axs[row,col].set_title(f'{case}')
    axs[row,col].set_xlabel('Density',fontsize=14)
    axs[row,col].set_ylabel('Difference in crosslink\n distances between native\n structure and models (Å)',fontsize=14)
    axs[row,col].tick_params(axis='both', which='major', labelsize=12)
    axs[row, col].legend(handles=[mpatches.Patch(color='blue'), mpatches.Patch(color='orange')], labels=['IMP', 'EASAL'])

plt.savefig(f'/home/muskaan/easal/plots/distance_distribution/{sys.argv[1]}.png')
plt.show()
