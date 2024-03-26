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
        #TODO very inefficient implementation.
        # Loading frame for each residue in the file
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

        #TODO storing heavy lists of objects that are not useful downstream. No need to store the core XYZ object if you are already storing one xlink per class
        #TODO

    def set_xl_coords(
        self, hierarchy: IMP.atom.Hierarchy, frame_id: int, rmf: RMF.FileConstHandle
    ) -> None:
        if len(self.xl.split(",")) != 4:
            raise ValueError(
                "Crosslinks file should have only four columns - r1, p1, r2, p2 - in this order"
            )

        r1, p1, r2, p2 = self.xl.split(",")

        #TODO this is not particle but a residue. Rename variable accordingly
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

        #TODO why do you need to call a function to load every frame and get the coords?
        #TODO can call it from right here. IMP.core.XYZ(prot1).get_coords()

        prot1.set_coords()
        prot2.set_coords()
        self.p1_xyz = prot1.coords
        self.p2_xyz = prot2.coords

    def get_min_distance(self) -> list:

    	#TODO these are binary complexes so we dont need min_distance here
        distances = []
        for c1 in self.p1_xyz:
            for c2 in self.p2_xyz:
                dist = IMP.core.get_distance(c1, c2)
                if dist < 0:
                    dist = 0
                distances.append(dist)

        return distances
        #TODO expected to return bool but returns a list of distances instead


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

xlink_filename = os.path.splitext(os.path.basename(xl_file))[0]
xl_dist = os.path.join('/home/muskaan/easal_imp/crosslink_distances/', f'{xlink_filename}_distances.txt')
percentage_satisfied = os.path.join('/home/muskaan/easal_imp/xl_satisfaction/', f'{xlink_filename}_percentage_satisfied.txt')

with open(xl_dist, 'w') as output_file, open(percentage_satisfied, 'w') as perc_satisfied:
    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        perc = 0
        for xl in all_xls:
            xl.set_xl_coords(hierarchy=hier, frame_id=frame, rmf=rmf_fh)
            xl_valid = xl.get_min_distance()
            output_file.write(f'{frame} {xl_valid[0]}\n')
            if xl_valid[0] < threshold:
                perc += 1

        #   #TODO return the actual xlink distance as we will need it in some of the plots. Even if violated.        #
            #TODO figure out a way to store xlink distances from each file such that the raw data from one script can be used for all the questions we are interested in.
            #TODO #framenum #crosslinkdist work ? DOnt add unnecessary words like frame and ':' . Can potentially use pandas or other tools to plot from these files.


        perc_satisfied.write(f'{frame} {perc/len(all_xls) *100}\n')
