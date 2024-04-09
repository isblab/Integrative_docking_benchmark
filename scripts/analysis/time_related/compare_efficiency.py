import os, sys
import matplotlib.pyplot as plt

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

fig, ax = plt.subplots(figsize=(8, 8))
# print(time_points_imp, time_points_easal)
# plt.hist(time_points_imp)
# plt.hist(time_points_easal)
# plt.xlabel('Runtime (minutes)', fontsize=14)
# plt.ylabel('Number of cases', fontsize=14)

plt.scatter(time_points_imp, time_points_easal)
plt.xlabel('Runtime for IMP (minutes)', fontsize=14)
plt.ylabel('Runtime for EASAL (minutes)', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
# plt.xlim(0, 3500)
# plt.ylim(0, 3500)
plt.savefig('/home/muskaan/easal/plots/time_taken_min.png',dpi=600)
plt.show()
