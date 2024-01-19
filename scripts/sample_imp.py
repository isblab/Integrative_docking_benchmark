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

def add_pdb_rep(mol,pdbname,chain,clr):
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

def getSequenceFromPDB(pdbfile,chain):
        seq = ""
        ms = IMP.Model()
        h = IMP.atom.read_pdb(pdbfile,ms)
        for ch in h.get_children():
             if ch!=chain:
                continue
	     for res in ch.get_children():
	         seq = seq + getAA_Alphabet(res.get_name())

	return seq

def recenter(mol):
    "recenter the system using the coordinate of a chain specified by its chain_id"
    m=mol.mdl
    sel = IMP.atom.Selection(mol.get_hierarchy(),resolution=IMP.atom.ALL_RESOLUTIONS)
    ps = sel.get_selected_particles()

    rb=IMP.atom.create_rigid_body(ps)
    rbcoord=rb.get_coordinates()

    rot=IMP.algebra.get_identity_rotation_3d()
    tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
    trans=tmptrans.get_inverse()
    IMP.core.transform(rb,trans)
    IMP.core.RigidBody.teardown_particle(rb)
    m.remove_particle(rb.get_particle_index())
    # change this back!Now rigid member coords need to be optimized as no longer in rigid body
    for p in ps:
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)

    return(mol)

def get_radius_center_residue(mol):

    sel = IMP.atom.Selection(mol.get_hierarchy(),resolution=1)
    ps = sel.get_selected_particles()

    max_dist=0.0
    min_dist=100000.0

    for p in ps:
        crd= IMP.core.XYZ(p).get_coordinates()

        dist=IMP.algebra.get_distance(crd,IMP.algebra.Vector3D(0.,0.,0.))

        if dist<min_dist:
            min_dist=dist
            center_res=p

        if dist>max_dist:
            max_dist=dist

    return(max_dist,center_res)

###################### SYSTEM SETUP #####################
# Parameters to tweak
tgt = sys.argv[1]
runType = sys.argv[2] # "prod" for production run and "test" for test run
inputDir=sys.argv[3]

xlink_file = inputDir+tgt+'.xlinks'

rpdb = inputDir+tgt+'_r_u.pdb'
rpdb_chain="A" #TODO command line args
lpdb = inputDir+tgt+'_l_u.pdb'
lpdb_chain="B"

receptor_chain = "A"
receptor_color = "blue"

ligand_chain = "B"
ligand_color = "red"


# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()


# Add Molecules for each component as well as representations
receptor = st.create_molecule(receptor_chain,chain_id=receptor_chain,sequence=getSequenceFromPDB(rpdb,rpdb_chain)) #TODO take IMP function
receptor = add_pdb_rep(receptor,rpdb,rpdb_chain,receptor_color)
ligand = st.create_molecule(ligand_chain,chain_id=ligand_chain,sequence=getSequenceFromPDB(lpdb,lpdb_chain))
ligand = add_pdb_rep(ligand,lpdb,lpdb_chain,ligand_color)


mols = [receptor,ligand]
###  Once you call build(), anything without representation is destroyed.
###  You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
###  However these functions will only return BUILT representations
root_hier = s.build()

RB_MAX_TRANS = 5.0 #TODO optimize
RB_MAX_ROT = 0.1

### Uncomment this for verbose output of the representation
#IMP.atom.show_with_representations(root_hier)


### output to RMF
#fname = 'testbefore.rmf'
#rh = RMF.create_rmf_file(fname)
#IMP.rmf.add_hierarchy(rh, root_hier)
#IMP.rmf.save_frame(rh)

receptor=recenter(receptor) #TODO may not need 
ligand=recenter(ligand)

(receptor_radius,receptor_center_residue)=get_radius_center_residue(receptor)
(ligand_radius,ligand_center_residue)=get_radius_center_residue(ligand)

print "radius",receptor_radius +ligand_radius+5.0

#fname = 'testafter.rmf'
#rh = RMF.create_rmf_file(fname)
#IMP.rmf.add_hierarchy(rh, root_hier)
#IMP.rmf.save_frame(rh)

## Setup degrees of freedom
##  The DOF functions automatically select all resolutions
##  Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
dof.create_rigid_body(ligand,max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,name="ligand")

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
evr.set_weight(0.03)
evr.add_to_model()
output_objects.append(evr)

# Crosslink restraint
# Not using the proxl database loader for now
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb.create_set_from_file(xlink_file)

#xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,CrossLinkDataBase=xldb,length=18.0,label="XLDSS",filelabel='dss',resolution=1,slope=0.05)
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,CrossLinkDataBase=xldb,length=8.0,label="XLDSS",filelabel='dss',resolution=1,slope=0.05)
xlr.add_to_model()
output_objects.append(xlr)
display_restraints.append(xlr)

# add distance to point restraint from receptor to ligand, to keep ligand floating away
#TODO try not to have this. See if ligand floats away.
ligand_center_residue_index=IMP.atom.Residue(ligand_center_residue).get_index()
receptor_center_residue_coords=IMP.core.XYZ(receptor_center_residue).get_coordinates()
dsr= IMP.pmi.restraints.basic.DistanceToPointRestraint(root_hier=root_hier,tuple_selection=(ligand_center_residue_index,ligand_center_residue_index,ligand_chain,0),anchor_point=receptor_center_residue_coords,radius=receptor_radius + ligand_radius + 3.0,weight=3.0)
dsr.add_to_model()
output_objects.append(dsr)

######################## SAMPLING #####################
## First shuffle the system
#IMP.pmi.tools.shuffle_configuration(root_hier.get_children()[0].get_children()[1],max_translation=500)

if runType == "prod":
    nframes = 30000

elif runType == "test":
    nframes = 600

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    crosslink_restraints=display_restraints,                     # will display like XLs
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    num_sample_rounds = 1,
                                    number_of_best_scoring_models = 500,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=nframes
				                    )

rex.execute_macro()
