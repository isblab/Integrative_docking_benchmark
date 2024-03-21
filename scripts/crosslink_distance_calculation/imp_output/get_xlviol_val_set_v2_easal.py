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


class Particle:
    """
    Constructor class for a particle. Contains the protein name, residue id and coordinates of corresponding beads
    """

    def __init__(
        self,
        prot: str,
        residue: int,
        hierarchy: IMP.atom.Hierarchy,
        frame_id: int,
        rmf: RMF.FileConstHandle,
    ) -> None:
        self.rmf_fh = rmf
        self.frame_id = frame_id
        self.hier = hierarchy
        self.pname = prot
        self.resid = residue
        self.coords = []

    def set_coords(self) -> None:
        IMP.rmf.load_frame(self.rmf_fh, self.frame_id)
        sel0 = IMP.atom.Selection(
            hierarchy=self.hier,
            molecule=self.pname,
            residue_index=self.resid,
            resolution=1,
        ).get_selected_particles()  # TODO crosscheck if all ambiguous beads are selected
        for bead in sel0:
            self.coords.append(IMP.core.XYZ(bead))


# TODO did we test that this script gives the same distances as in stat file or RMF for ambiguous crosslinks?


class Xlink:
    """
    Constructor for a crosslink class. Reads a file and generates crosslinks with particles from the Particle class
    """

    def __init__(
        self,
        xl_ln: str,
    ) -> None:
        self.xl = xl_ln.strip()
        self.p1_xyz = []
        self.p2_xyz = []
        self.min_distances = []  # minimum distance from each model
        self.violated = None

    def set_xl_coords(
        self, hierarchy: IMP.atom.Hierarchy, frame_id: int, rmf: RMF.FileConstHandle
    ) -> None:
        if len(self.xl.split(",")) != 4:
            raise ValueError(
                "Crosslinks file should have only four columns - r1, p1, r2, p2 - in this order"
            )

        r1, p1, r2, p2 = self.xl.split(",")
        prot1 = Particle(
            prot=p1,
            residue=int(r1),
            hierarchy=hierarchy,
            frame_id=frame_id,
            rmf=rmf,
        )

        prot2 = Particle(
            prot=p2,
            residue=int(r2),
            hierarchy=hierarchy,
            frame_id=frame_id,
            rmf=rmf,
        )

        prot1.set_coords()
        prot2.set_coords()
        self.p1_xyz = prot1.coords
        self.p2_xyz = prot2.coords

    def get_min_distance(self) -> bool:
        distances = []
        for c1 in self.p1_xyz:
            for c2 in self.p2_xyz:
                dist = IMP.core.get_distance(c1, c2)
                if dist < 0:
                    dist = 0
                distances.append(dist)

        if len(distances) > 0:
            self.min_distances.append(min(distances))
            return distances
        # else:
        #     return False

    def set_violation_status(self, viol_threshold):
        self.violated = True
        for dist in self.min_distances:
            if dist < viol_threshold:
                self.violated = False

all_xls: list[Xlink] = []
with open(xl_file, "r") as xlf:
    for ln in xlf.readlines():
        if not ln.startswith("res1"):
            xl = Xlink(ln.strip())
            all_xls.append(xl)

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(all_mdls_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
mdl.update()

all_valid_xls: list[Xlink] = []
xlink_filename = os.path.splitext(os.path.basename(xl_file))[0]
xl_dist = os.path.join(os.getcwd(), f'{xlink_filename}_distances.txt')

with open(xl_dist, 'w') as output_file:
    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        frame_dist = []
        for xl in all_xls:
            xl.set_xl_coords(hierarchy=hier, frame_id=frame, rmf=rmf_fh)
            xl_valid = xl.get_min_distance()
            frame_dist.append((frame, xl_valid[0]) if xl_valid else (frame, None))
            if frame == 0 and xl_valid:
                all_valid_xls.append(xl)

        for frame_num, xl_valid_value in frame_dist:
            output_file.write(f'Frame {frame_num}: {xl_valid_value}\n')
