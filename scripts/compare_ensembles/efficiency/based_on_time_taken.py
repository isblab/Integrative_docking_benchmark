import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime


# easal

input_cases = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_DMTMM_4", "1clv_DMTMM_8", "1kxp_DMTMM_7", "1r0r_DMTMM_6", "2ayo_DMTMM_5", "2b42_DMTMM_10", "2hle_DMTMM_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]

time_points_easal, time_per_run, time_points_imp = [],[],[]

for case in input_cases:
    if 'DSSO' in case and len(case) <15:
        file = '/home/muskaan/easal/time_related/DSSO/simulated/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'

    elif 'DSSO' in case and len(case) >15:
        file = '/home/muskaan/easal/time_related/DSSO/experimental/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
    else:
        file = '/home/muskaan/easal/time_related/DMTMM/'+case.split('DMTMM')[0] + case.split('_')[-1]+ '/logfile.txt'

    with open(file, 'r') as logfile:
        for line in logfile.readlines():
            if 'Thread_Main: Total time cost for sampling' in line:
                time = (int(line.split('sampling:')[-1])/1000)/60 #In minutes (milisec/1000)/60
                time_points_easal.append(time)
                # print(case, time)
                # print(line.split('sampling:')[-1])


# imp
for case in input_cases:
    if 'DSSO' in case and len(case) <15:
        file1 = '/home/muskaan/easal/imp_output/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/run_1/stat_replica.0.out'
        file2 = '/home/muskaan/easal/imp_output/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/run_1/initial.0.rmf3'

    elif 'DSSO' in case and len(case) >15:
        file1 = '/home/muskaan/easal/imp_output/DSSO/'+case.split('_DSSO')[0] + '/run_1/stat_replica.0.out'
        file2 = '/home/muskaan/easal/imp_output/DSSO/'+case.split('_DSSO')[0] + '/run_1/initial.0.rmf3'

    elif 'DMTMM' in case:
        file1 = '/home/muskaan/easal/imp_output/DMTMM/'+case.split('DMTMM')[0] + case.split('_')[-1]+ '/run_1/stat_replica.0.out'
        file2 = '/home/muskaan/easal/imp_output/DMTMM/'+case.split('DMTMM')[0] + case.split('_')[-1]+ '/run_1/initial.0.rmf3'

    time = os.path.getmtime(file1) - os.path.getmtime(file2)
    time_total = (time/60) * 20 * 4 #In minutes; per run for 4 replica and 20 runs, multiply by 20 *4
    time_per_run.append(time/60)
    time_points_imp.append(time_total)

# print(np.mean(time_points_imp))
# print(np.mean(time_points_easal))
# print(np.mean(time_per_run))

methods = ['IMP', 'EASAL']
x_positions = np.arange(len(methods))+1

fig, ax = plt.subplots(figsize=(10, 8))
ax.bar(1, np.mean(time_points_imp), yerr=[[0], [np.std(time_points_imp)]], label = 'IMP', color ='#1f77b4', capsize=5)
ax.bar(2, np.mean(time_points_easal), yerr=[[0],[np.std(time_points_easal)]], label = 'EASAL', color ='#ff7f0e',capsize=5)
ax.set_xlabel('Method', fontsize=18)
ax.set_ylabel('Total runtime (CPU hours)', fontsize=18)
ax.set_xticks(x_positions)
ax.set_xticklabels(methods, fontsize=14)
plt.savefig('/home/muskaan/easal/plots/time_related/F6.runtime.png',dpi=600)
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
ax.bar(1, np.mean(time_per_run), yerr=[[0],[np.std(time_per_run)]], label = 'IMP', color ='#1f77b4',capsize=5)
ax.bar(2, np.mean(time_points_easal), yerr=[[0],[np.std(time_points_easal)]], label = 'EASAL', color ='#ff7f0e',capsize=5)
ax.set_xlabel('Method', fontsize=18)
plt.ylabel('Sampling time (minutes per run)', fontsize=18)
ax.set_xticks(x_positions)
ax.set_xticklabels(methods, fontsize=14)
plt.savefig('/home/muskaan/easal/plots/time_related/F6.runtime_per_run.png',dpi=600)
plt.show()
