import os
import numpy as np
import matplotlib.pyplot as plt

def get_imp_easal_efficiency_data(input_cases, fp=True):
    imp_ratio, easal_ratio, plot_cases = [],[],[]
    total_easal, total_imp = [], []
    # IMP
    for case in input_cases:
        if fp:
            xl_sat_file_imp = '/home/muskaan/projects/easal_related/easal/imp_output/test_fp/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
        else:
            xl_sat_file_imp = '/home/muskaan/projects/easal_related/easal/imp_output/xl_satisfaction/' + f'{case}_perc_satisfied.txt'

        total_perc_imp = []
        with open(xl_sat_file_imp, 'r') as xl_sat_file_imp:
            for line in xl_sat_file_imp:
                models_perc = float(line.split()[1])
                total_perc_imp.append(models_perc)

        max_xlink_perc_imp = max(total_perc_imp)
        models_with_max_xlink_sat_imp = sum(1 for perc in total_perc_imp if perc >= max_xlink_perc_imp)
        # print(max_xlink_perc_imp, case)

        ## EASAL
        if fp:
            xl_sat_file_easal = '/home/muskaan/projects/easal_related/easal/easal_output/with_fp/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
            file = '/home/muskaan/projects/easal_related/easal/easal_output/with_fp/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
        else:
            # Using newer EASAL results
            xl_sat_file_easal = '/home/muskaan/projects/easal_related/easal/easal_output/without_fp/xl_satisfaction/' + f'{case}_perc_satisfied.txt'
            if 'DSSO' in case:
                file = '/home/muskaan/projects/easal_related/easal/easal_output/without_fp/DSSO/'+case.split('DSSO')[0] + case.split('_')[-1]+ '/logfile.txt'
            elif 'DMTMM' in case:
                file = '/home/muskaan/projects/easal_related/easal/easal_output/without_fp/DMTMM/'+case.split('DMTMM')[0] + case.split('_')[-1]+ '/logfile.txt'

            # Using older EASAL results
            xl_sat_file_easal = '/home/muskaan/projects/easal_related/easal/easal_output/older_version/xl_satisfaction/' + f'{case}_perc_satisfied.txt'

        total_perc_easal = []

        with open(xl_sat_file_easal, 'r') as xl_sat_file_easal:
            for line in xl_sat_file_easal:
                models_perc = float(line.split()[1])
                total_perc_easal.append(models_perc)

        max_xlink_perc_easal = max(total_perc_easal)

        if max_xlink_perc_imp == max_xlink_perc_easal: # Plot if the max xlink sat is same of IMP and wall-EASAL ensembles
            with open(file, 'r') as logfile:
                lines = logfile.readlines()
                for i, line in enumerate(lines):
                    if 'Best sample count:' in line:
                        total_sample_easal = int(lines[i - 1].split('count:')[-1])
                        best_sample_count = int(lines[i].split('count:')[-1])
                        # print(case , total_sample_easal, best_sample_count)

            easal_ratio.append(best_sample_count/total_sample_easal)
            imp_ratio.append(models_with_max_xlink_sat_imp/8000000)
            plot_cases.append(case)
            total_easal.append(total_sample_easal)
            total_imp.append(8000000)
        # print(models_with_max_xlink_sat_easal, easal_ratio, case)
        # print(models_with_max_xlink_sat_imp, imp_ratio, case)

    # easal_ratio = [ratio * 100 for ratio in easal_ratio]
    # imp_ratio = [ratio * 100 for ratio in imp_ratio]

    return plot_cases, imp_ratio, easal_ratio, total_easal, total_imp

input_cases_fp = [ "1clv_DSSO_5", "1dfj_DSSO_10", "1dfj_DSSO_11", "1kxp_DSSO_10", "1r0r_DSSO_5", "2ayo_DSSO_10", "2ayo_DSSO_15", "2b42_DSSO_10", "2hle_DSSO_10", "2hle_DSSO_15"]

input_cases_tp = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
    "1dfj_DMTMM_4", "1clv_DMTMM_8", "1kxp_DMTMM_7", "1r0r_DMTMM_6", "2ayo_DMTMM_5", "2b42_DMTMM_10", "2hle_DMTMM_9",
    "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
    "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
    "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]


plot_cases_fp, imp_fp, easal_fp, total_easal_fp, total_imp_fp = get_imp_easal_efficiency_data(input_cases_fp, True)
plot_cases_tp, imp_tp, easal_tp, total_easal_tp, total_imp_tp = get_imp_easal_efficiency_data(input_cases_tp, False)

# plt.figure(figsize = (40,40))
fig, axs = plt.subplots(2, 2, figsize = (15,15))

axs[0, 0].scatter(plot_cases_tp, easal_tp, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
axs[0, 0].scatter(plot_cases_tp, imp_tp, color='#1f77b4', label='IMP', alpha=0.7)
axs[0, 0].set_xlabel('Input Cases (without false positives)',  fontsize=14)
axs[0, 0].set_ylabel('Fraction of best\n configurations in the sample',  fontsize=14)
axs[0, 0].set_xticks(plot_cases_tp)
axs[0, 0].set_xticklabels(plot_cases_tp, rotation=45, ha='right', fontsize=8)

axs[0, 0].legend()

axs[0, 1].scatter(plot_cases_tp, total_easal_tp, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
axs[0, 1].scatter(plot_cases_tp, total_imp_tp, color='#1f77b4', label='IMP', alpha=0.7)
axs[0, 1].set_xlabel('Input Cases (without false positives)',  fontsize=14)
axs[0, 1].set_ylabel('Total samples',  fontsize=14)
axs[0, 1].set_xticks(plot_cases_tp)
axs[0, 1].set_xticklabels(plot_cases_tp, rotation=45, ha='right', fontsize=8)
axs[0, 1].legend()

axs[1, 0].scatter(plot_cases_fp, easal_fp, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
axs[1, 0].scatter(plot_cases_fp, imp_fp, color='#1f77b4', label='IMP', alpha=0.7)
axs[1, 0].set_xlabel('Input Cases (with false positives)',  fontsize=14)
axs[1, 0].set_ylabel('Fraction of best\n configurations in the sample',  fontsize=14)
axs[1, 0].set_xticks(plot_cases_fp)
axs[1, 0].set_xticklabels(plot_cases_fp, rotation=45, ha='right', fontsize=8)
axs[1, 0].legend()


axs[1, 1].scatter(plot_cases_fp, total_easal_fp, color='#ff7f0e', label='Wall-EASAL', alpha=0.7)
axs[1, 1].scatter(plot_cases_fp, total_imp_fp, color='#1f77b4', label='IMP', alpha=0.7)
axs[1, 1].set_xlabel('Input Cases (with false positives)',  fontsize=14)
axs[1, 1].set_ylabel('Total samples',  fontsize=14)
axs[1, 1].set_xticks(plot_cases_fp)
axs[1, 1].set_xticklabels(plot_cases_fp, rotation=45, ha='right', fontsize=8)
axs[1, 1].legend()

fig.subplots_adjust(hspace=0.4, wspace=0.4)
plt.tight_layout()
plt.savefig('/home/muskaan/projects/easal_related/easal/plots/time_related/F10.sampling_efficiency_newer.png',dpi=600)
plt.show()
