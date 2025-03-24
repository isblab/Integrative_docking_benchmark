import os,sys,math
import numpy as np
import glob
import csv
import matplotlib.pyplot as plt

input_cases = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5","2hle_DSSO_5",
    "1dfj_DMTMM_4", "1clv_DMTMM_8", "1kxp_DMTMM_7", "1r0r_DMTMM_6", "2ayo_DMTMM_5",
    "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10", "1dfj_DSSO_9",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "sucd_succ_DSSO_4", "roca_putc_DSSO_2"]

## IMP
imp_ratio, easal_ratio = [],[]

for case in input_cases:
    xl_sat_file = '~/easal/imp_output/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
    total_perc = []
    with open(xl_sat_file, 'r') as xl_sat_file:
        for line in xl_sat_file:
            models_perc = float(line.split()[1])
            total_perc.append(models_perc)

    max_xlink_perc = max(total_perc)
    models_with_max_xlink_sat = sum(1 for perc in total_perc if perc >= max_xlink_perc)

    imp_ratio.append(models_with_max_xlink_sat/8000000)
    # print(models_with_max_xlink_sat, imp_ratio, case)

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
    # print(file_count/total_sample_count, case)
    easal_ratio.append(file_count/total_sample_count)

easal_ratio = [ratio * 100 for ratio in easal_ratio]
imp_ratio = [ratio * 100 for ratio in imp_ratio]

print(imp_ratio, easal_ratio )
plt.figure(figsize=(12, 6))
plt.scatter(input_cases, easal_ratio, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
plt.scatter(input_cases, imp_ratio, color='#1f77b4', label='IMP', alpha=0.7)

plt.xlabel('Input Cases',  fontsize=14)
plt.ylabel('Fraction of best configurations in the sample',  fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.text(0.01, 1.05, r'$x10^{-2}$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.tight_layout()
plt.savefig('~/easal/plots/time_related/F10.sampling_efficiency.png',dpi=600)
plt.show()
