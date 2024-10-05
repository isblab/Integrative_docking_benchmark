import sys, os
import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom
from tqdm import tqdm

all_mdls_fname = sys.argv[1]
xl_file = sys.argv[2]
threshold = float(sys.argv[3])


def calculate_distance(xl):
    r1, p1, r2, p2 = xl.split(',')
    c1 = IMP.atom.Selection(hierarchy=hier,
                            molecule=p1,
                            resolution=1,
                            residue_index=int(r1)).get_selected_particles()

    c2 = IMP.atom.Selection(hierarchy=hier,
                            molecule=p2,
                            resolution=1,
                            residue_index=int(r2)).get_selected_particles()

    dist = IMP.core.get_distance(IMP.core.XYZ(c1[0]), IMP.core.XYZ(c2[0]))
    return dist

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(all_mdls_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
all_xls = []

with open(xl_file, "r") as xlf:
    for ln in xlf.readlines():
        if not ln.startswith("res1"):
            xl = ln.strip()
            all_xls.append(xl)


xlink_filename = os.path.splitext(os.path.basename(xl_file))[0]
xl_dist = os.path.join('easal/imp_output/crosslink_distances/', f'{xlink_filename}_distances.txt')
percentage_satisfied = os.path.join('easal/imp_output/xl_satisfaction/', f'{xlink_filename}_perc_satisfied.txt')

with open(xl_dist, 'w') as output_dist_file, open(percentage_satisfied, 'w') as output_perc_sat_file:

    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_fh, frame)
        mdl.update()
        perc = 0
        for xl in all_xls:
            distance = calculate_distance(xl)
            if distance < 0:
                distance = 0
            output_dist_file.write(f'{frame} {distance}\n')
            if distance < threshold:
                perc += 1
        output_perc_sat_file.write(f'{frame} {(perc/len(all_xls)) *100}\n')
