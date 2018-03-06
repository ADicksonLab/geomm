import numpy as np
import mdtraj as mdj

from geomm.recentering import group_pair, apply_rectangular_pbcs

# helper functions for getting ligand and protein idxs
def ligand_idxs(mdtraj_topology, ligand_resid):
    return mdtraj_topology.select('resname "{}"'.format(ligand_resid))

def protein_idxs(mdtraj_topology):
    return np.array([atom.index for atom in mdtraj_topology.atoms if atom.residue.is_protein])

uncentered_system = mdj.load_pdb("tryp_ben_system.pdb")
top = uncentered_system.topology

lig_idxs = ligand_idxs(top, 'BEN')
prot_idxs = protein_idxs(top)

# This particular frame of coordinates comes from an OpenMM simulation
# which has it's box centered at the origin
box_center = np.array([0., 0., 0.,])


# one thing we can do is wrap atoms around into the box
clipped_coords = apply_rectangular_pbcs(uncentered_system.xyz[0],
                                        uncentered_system.unitcell_lengths[0],
                                        box_center)
# write this out to check it
clipped_system = mdj.Trajectory(clipped_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

clipped_system.save_pdb("clipped_tryp_ben.pdb")


# first group the protein and ligand together, this may be outside of
# the PBCs
grouped_coords = group_pair(uncentered_system.xyz[0], uncentered_system.unitcell_lengths[0],
                            prot_idxs, lig_idxs)

# write this out to check it
grouped_system = mdj.Trajectory(grouped_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

grouped_system.save_pdb("grouped_tryp_ben.pdb")



# now we want to recenter the system around the protein-ligand complex
# so that it can be in the center and not have to be outside the box
recentered_coords = (grouped_coords,
                                           uncentered_system.unitcell_lengths[0],
                                           box_center)

# write them out so we can visualize
recentered_system = mdj.Trajectory(recentered_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

recentered_system.save_pdb("recentered_tryp_ben.pdb")

