import os
import numpy as np
import matplotlib.pyplot as plt

input_cases = [ "1clv_DSSO_5", "1dfj_DSSO_10", "1dfj_DSSO_11", "1kxp_DSSO_10", "1r0r_DSSO_5", "2ayo_DSSO_10", "2ayo_DSSO_15", "2b42_DSSO_10", "2hle_DSSO_10", "2hle_DSSO_15"]

## IMP
imp_ratio, easal_ratio, plot_cases = [],[],[]

for case in input_cases:
    xl_sat_file_imp = '~/easal/imp_output/test_fp/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
    total_perc_imp = []
    with open(xl_sat_file_imp, 'r') as xl_sat_file_imp:
        for line in xl_sat_file_imp:
            models_perc = float(line.split()[1])
            total_perc_imp.append(models_perc)

    max_xlink_perc_imp = max(total_perc_imp)
    models_with_max_xlink_sat_imp = sum(1 for perc in total_perc_imp if perc >= max_xlink_perc_imp)
    # print(max_xlink_perc_imp, case)

    ## EASAL
    # for case in input_cases:
    xl_sat_file_easal = '~/easal/easal_output/test_fp/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
    file = '~/easal/easal_output/test_fp/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'

    with open(file, 'r') as logfile:
        lines = logfile.readlines()
        for i, line in enumerate(lines):
            if 'FINISHED AtlasBuilder::startAtlasBuilding' in line:
                total_sample_count = int(lines[i - 1].split('count:')[-1])
                # print(case , total_sample_count)

    total_perc_easal = []
    with open(xl_sat_file_easal, 'r') as xl_sat_file_easal:
        for line in xl_sat_file_easal:
            models_perc = float(line.split()[1])
            total_perc_easal.append(models_perc)

    max_xlink_perc_easal = max(total_perc_easal)

    if max_xlink_perc_imp == max_xlink_perc_easal: # Plot if the max xlink sat is same of IMP and wall-EASAL ensembles
        models_with_max_xlink_sat_easal = sum(1 for perc in total_perc_easal if perc >= max_xlink_perc_easal)
        # print(max_xlink_perc_easal, case)
        easal_ratio.append(models_with_max_xlink_sat_easal/total_sample_count)
        imp_ratio.append(models_with_max_xlink_sat_imp/8000000)
        plot_cases.append(case)
    # print(models_with_max_xlink_sat_easal, easal_ratio, case)
    # print(models_with_max_xlink_sat_imp, imp_ratio, case)

easal_ratio = [ratio * 10000 for ratio in easal_ratio]
imp_ratio = [ratio * 10000 for ratio in imp_ratio]

print(imp_ratio, easal_ratio)
plt.figure(figsize=(12, 6))
plt.scatter(plot_cases, easal_ratio, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
plt.scatter(plot_cases, imp_ratio, color='#1f77b4', label='IMP', alpha=0.7)

plt.xlabel('Input Cases',  fontsize=14)
plt.ylabel('Fraction of best configurations in the sample',  fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.text(0.01, 1.05, r'$x10^{-4}$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.tight_layout()
plt.savefig('~/easal/plots/time_related/F10.sampling_efficiency_fp.png',dpi=600)
plt.show()
