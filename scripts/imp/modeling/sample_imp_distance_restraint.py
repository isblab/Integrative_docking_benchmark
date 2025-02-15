import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import os,sys,math,numpy
import ihm
import glob

def add_pdb_rep(mol,pdbname,chain,clr):
    print(chain)
    atomic = mol.add_structure(pdbname,chain_id=chain,soft_check=True)
    mol.add_representation(atomic, resolutions=[1],color = clr)
    return mol

def getAA_Alphabet(resnam):

	if resnam=='ALA':
		return 'A'
	elif resnam=='CYS':
		return 'C'
	elif resnam=='ASP':
		return 'D'
	elif resnam=='GLU':
		return 'E'
	elif resnam=='PHE':
		return 'F'
	elif resnam=='GLY':
		return 'G'
	elif resnam=='HIS':
		return 'H'
	elif resnam=='ILE':
		return 'I'
	elif resnam=='LYS':
		return 'K'
	elif resnam=='LEU':
		return 'L'
	elif resnam=='MET':
		return 'M'
	elif resnam=='ASN':
		return 'N'
	elif resnam=='PRO':
		return 'P'
	elif resnam=='GLN':
		return 'Q'
	elif resnam=='ARG':
		return 'R'
	elif resnam=='SER':
		return 'S'
	elif resnam=='THR':
		return 'T'
	elif resnam=='VAL':
		return 'V'
	elif resnam=='TRP':
		return 'W'
	elif resnam=='TYR':
		return 'Y'
	else:
		return 'G'
	# default is GLYcine

def getSequenceFromPDB(pdbfile, chain):
    seq = ""
    ms = IMP.Model()
    h = IMP.atom.read_pdb(pdbfile, ms)
    for ch in h.get_children():
        if ch != chain:
            continue
        for res in ch.get_children():
            seq += getAA_Alphabet(res.get_name())
    return seq

###################### SYSTEM SETUP #####################
# Parameters to tweak

runType = sys.argv[1] # "prod" for production run and "test" for test run
rpdb_chain = sys.argv[2]
lpdb_chain = sys.argv[3]
run_output_dir = 'run_' + sys.argv[4]
crosslinker = sys.argv[5]
xl_length = float(sys.argv[6])
xlink_file = glob.glob('*.csv')[0]

pdb = glob.glob('*.pdb')[0]

receptor_color = "blue"
ligand_color = "red"


# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

# Add Molecules for each component as well as representations
receptor = st.create_molecule(rpdb_chain,chain_id=rpdb_chain,sequence=getSequenceFromPDB(pdb,rpdb_chain)) #TODO take IMP function
receptor = add_pdb_rep(receptor,pdb,rpdb_chain,receptor_color)
ligand = st.create_molecule(lpdb_chain,chain_id=lpdb_chain,sequence=getSequenceFromPDB(pdb,lpdb_chain))
ligand = add_pdb_rep(ligand,pdb,lpdb_chain,ligand_color)

mols = [receptor,ligand]
###  Once you call build(), anything without representation is destroyed.
###  You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
###  However these functions will only return BUILT representations
root_hier = s.build()

RB_MAX_TRANS = 5.0
RB_MAX_ROT = 0.5

## Setup degrees of freedom
##  The DOF functions automatically select all resolutions
##  Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
dof.create_rigid_body(ligand,max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,name="ligand")
dof.create_rigid_body(receptor,max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,name="receptor")

######################## RESTRAINTS #####################
output_objects = [] # keep a list of functions that need to be reported
display_restraints = [] # display as springs in RMF

# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
for mol in mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    output_objects.append(cr)
    crs.append(cr)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.set_weight(0.5)
evr.add_to_model()
output_objects.append(evr)

# Harmonic distance restraint for crosslink data
xl = []
with open(xlink_file, 'r') as xlfile:
    for line in xlfile.readlines()[1:]:
        xl.append(((int(line.strip().split(',')[0]), int(line.strip().split(',')[0]), str(rpdb_chain),0), (int(line.strip().split(',')[2]), int(line.strip().split(',')[2]), str(lpdb_chain),0)))

for t in xl:
    xlr = IMP.pmi.restraints.basic.DistanceRestraint(root_hier=root_hier, tuple_selection1=t[0], tuple_selection2=t[1], distancemax=xl_length, distancemin=10)
    xlr.add_to_model()

output_objects.append(xlr)
display_restraints.append(xlr)

######################## SAMPLING #####################
## First shuffle the system
IMP.pmi.tools.shuffle_configuration(root_hier.get_children()[0].get_children()[1],max_translation=50)

if runType == "prod":
    nframes = 10000

elif runType == "test":
    nframes = 100

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    # crosslink_restraints=display_restraints,                     # will display like XLs
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory=run_output_dir,
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=nframes
				                    )

rex.execute_macro()
