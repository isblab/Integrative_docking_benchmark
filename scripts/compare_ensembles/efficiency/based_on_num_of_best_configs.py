import os,sys,math
import numpy as np
import glob
import csv
import matplotlib.pyplot as plt

input_cases = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5",
    "1dfj_DMTMM_4", "1clv_DMTMM_8", "1kxp_DMTMM_7", "1r0r_DMTMM_6", "2ayo_DMTMM_5", "2hle_DMTMM_9",
    "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]


## IMP
imp_ratio, easal_ratio = [],[]
for case in input_cases:
    if 'DSSO' in case:
        base_path = '~/easal/imp_output/DSSO_analysis/' + case.split('DSSO')[0] + case.split('_')[-1] + '/model_analysis/'
        sample_A_file = base_path + 'A_models*.txt'
        sample_B_file = base_path + 'B_models*.txt'
    else:
        base_path = '~/easal/imp_output/DMTMM_analysis/' + case.split('DMTMM')[0] + case.split('_')[-1] + '/model_analysis/'
        sample_A_file = base_path + 'A_models*.txt'
        sample_B_file = base_path + 'B_models*.txt'

    A_files = glob.glob(sample_A_file)
    num_rows_A = sum(sum(1 for line in open(a_file, 'r')) for a_file in A_files)

    B_files = glob.glob(sample_B_file)
    num_rows_B = sum(sum(1 for line in open(b_file, 'r')) for b_file in B_files)

    # print(num_rows_A + num_rows_B)

    xl_sat_file = '~/easal/imp_output/xl_satisfaction/' + f'{case}_perc_satisfied.txt'

    total_perc = []
    with open(xl_sat_file, 'r') as xl_sat_file:
        for line in xl_sat_file:
            models_perc = float(line.split()[1])
            total_perc.append(models_perc)

    max_xlink_perc = max(total_perc)
    models_with_max_xlink_sat = sum(1 for perc in total_perc if perc >= max_xlink_perc)

    # print(models_with_max_xlink_sat)
    print(models_with_max_xlink_sat/(num_rows_A + num_rows_B), case)
    imp_ratio.append(models_with_max_xlink_sat/(num_rows_A + num_rows_B))
    # imp_ratio.append(models_with_max_xlink_sat/8000000)

## EASAL
for case in input_cases:
    if 'DSSO' in case and len(case) <15:
        file = '~/easal/time_related/DSSO/simulated/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
        directory_path = '~/easal/easal_output/DSSO/simulated/' + case.split('DSSO')[0] + 'cl' + case.split('_')[-1]

    elif 'DSSO' in case and len(case) >15:
        file = '~/easal/time_related/DSSO/experimental/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
        directory_path = '~/easal/easal_output/DSSO/experimental/' + case.split('DSSO')[0] + 'cl' + case.split('_')[-1]
    else:
        file = '~/easal/time_related/DMTMM/'+case.split('DMTMM')[0] + case.split('_')[-1]+ '/logfile.txt'
        directory_path = '~/easal/easal_output/DMTMM/' + case.split('DMTMM')[0] + 'cl' + case.split('_')[-1]

    with open(file, 'r') as logfile:
        lines = logfile.readlines()
        for i, line in enumerate(lines):
            if 'FINISHED AtlasBuilder::startAtlasBuilding' in line:
                total_sample_count = int(lines[i - 1].split('count:')[-1])
                # print(case , total_sample_count)


    file_count = sum(1 for entry in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, entry)))
    # print(file_count, case)
    print(file_count/total_sample_count, case)
    easal_ratio.append(file_count/total_sample_count)



methods = ['IMP', 'EASAL']
x_positions = np.arange(len(methods))+1

fig, ax = plt.subplots(figsize=(10, 8))
ax.bar(1, np.mean(imp_ratio), yerr=[[0], [np.std(imp_ratio)]], label = 'IMP', color ='#1f77b4', capsize=5)
ax.bar(2, np.mean(easal_ratio), yerr=[[0],[np.std(easal_ratio)]], label = 'EASAL', color ='#ff7f0e',capsize=5)
ax.set_xlabel('Method', fontsize=18)
ax.set_ylabel('Fraction of best configurations in the sample', fontsize=18)
ax.set_xticks(x_positions)
ax.set_xticklabels(methods, fontsize=14)
plt.savefig('~/easal/plots/time_related/F6.efficiency.png',dpi=600)
# plt.savefig('~/easal/plots/time_related/F6.efficiency_total.png',dpi=600)
plt.show()
