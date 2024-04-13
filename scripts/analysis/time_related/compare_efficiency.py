import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# easal

input_cases = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]

time_points_easal, time_points_imp = [],[]

for case in input_cases:
    if 'DSSO' in case and len(case) <15:
        file = '/home/muskaan/easal/time_related/DSSO/simulated/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'

    elif 'DSSO' in case and len(case) >15:
        file = '/home/muskaan/easal/time_related/DSSO/experimental/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
    else:
        file = '/home/muskaan/easal/time_related/EDC/'+case.split('EDC')[0] + case.split('_')[-1]+ '/logfile.txt'

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
        file = '/home/muskaan/easal/imp_output/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/run_1/stat_replica.0.out'

    elif 'DSSO' in case and len(case) >15:
        file = '/home/muskaan/easal/imp_output/DSSO/'+case.split('_DSSO')[0] + '/run_1/stat_replica.0.out'

    elif 'EDC' in case:
        file = '/home/muskaan/easal/imp_output/EDC/'+case.split('EDC')[0] + case.split('_')[-1]+ '/run_1/stat_replica.0.out'

    time = os.path.getmtime(file)-os.path.getatime(file)
    time_total = (time/60) * 20 * 4 #In minutes; per run for 4 replica and 20 runs, multiply by 20 *4
    time_points_imp.append(time_total)

    # print(case, time_total)
fig, ax = plt.subplots(figsize=(10, 8))
plt.violinplot(time_points_imp, showmeans=False, showmedians=False)
plt.violinplot(time_points_easal, showmeans=False, showmedians=False)
plt.xlabel('Density', fontsize=18)
plt.ylabel('Runtime (minutes)', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(handles=[mpatches.Patch(color='#1f77b4'), mpatches.Patch(color='#ff7f0e')], labels=['IMP', 'EASAL'])
plt.savefig('/home/muskaan/easal/plots/time_related/F6.runtime.png',dpi=600)
plt.show()
