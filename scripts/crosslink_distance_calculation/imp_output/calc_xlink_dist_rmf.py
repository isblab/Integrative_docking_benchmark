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


def set_coords(xl):
    r1, p1, r2, p2 = xl.split(',')
    c1 = IMP.atom.Selection(hierarchy=hier,
                            molecule=p1,
                            resolution=1,
                            residue_index=int(r1)).get_selected_particles()

    c2 = IMP.atom.Selection(hierarchy=hier,
                            molecule=p2,
                            resolution=1,
                            residue_index=int(r2)).get_selected_particles()


    return c1[0], c2[0]

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
xl_dist = os.path.join('/home/muskaan/easal_imp/crosslink_distances/', f'{xlink_filename}_distances.txt')
percentage_satisfied = os.path.join('/home/muskaan/easal_imp/xl_satisfaction/', f'{xlink_filename}_percentage_satisfied.txt')

with open(xl_dist, 'w') as output_file, open(percentage_satisfied, 'w') as perc_satisfied:
    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_fh, frame)
        mdl.update()
        perc = 0
        for xl in all_xls:
            c1, c2 = set_coords(xl)
            distance = IMP.core.get_distance(IMP.core.XYZ(c1), IMP.core.XYZ(mdl, c2))
            if distance <0:
                distance = 0
            output_file.write(f'{frame} {distance}\n')
            if distance < threshold:
                perc += 1
        perc_satisfied.write(f'{frame} {perc/len(all_xls) *100}\n')
