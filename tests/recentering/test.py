import numpy as np
import mdtraj as mdj

from geomm.recentering import group_pair, apply_rectangular_pbcs, recenter_pair

# helper functions for getting ligand and protein idxs
def ligand_idxs(mdtraj_topology, ligand_resid):
    return mdtraj_topology.select('resname "{}"'.format(ligand_resid))

def protein_idxs(mdtraj_topology):
    return np.array([atom.index for atom in mdtraj_topology.atoms if atom.residue.is_protein])

uncentered_system = mdj.load_pdb("tryp_ben_system.pdb")
top = uncentered_system.topology
box_lengths = uncentered_system.unitcell_lengths[0]

lig_idxs = ligand_idxs(top, 'BEN')
prot_idxs = protein_idxs(top)

# This particular frame has its center at the (0.5*L_x, 0.5*L_y,
# 0.5*L_z), where the 'L's are the components of the unitcell lengths
box_center = np.array([box_lengths[i] * 0.5 for i in range(3)])

# this should wrap things into the box assuming it is centered as
# defined above
clipped_coords = apply_rectangular_pbcs(uncentered_system.xyz[0],
                                        box_lengths,
                                        box_center)
# write this out to check it
clipped_system = mdj.Trajectory(clipped_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

clipped_system.save_pdb("clipped_center_tryp_ben.pdb")



unit_corners = [(0,0,0), (0,0,1), (0,1,0), (1,0,0), (0,1,1), (1,0,1), (1,1,0), (1,1,1)]
for unit_corner in unit_corners:
    corner = unit_corner * box_lengths
    # one thing we can do is wrap atoms around into the box
    clipped_coords = apply_rectangular_pbcs(uncentered_system.xyz[0],
                                            box_lengths,
                                            corner)
    # write this out to check it
    clipped_system = mdj.Trajectory(clipped_coords, top,
                                    unitcell_lengths=uncentered_system.unitcell_lengths,
                                    unitcell_angles=uncentered_system.unitcell_angles)

    clipped_system.save_pdb("clipped_{}_{}_{}_tryp_ben.pdb".format(*unit_corner))


# first group the protein and ligand together, this may be outside of
# the PBCs
grouped_coords = group_pair(uncentered_system.xyz[0], box_lengths,
                            prot_idxs, lig_idxs)

# write this out to check it
grouped_system = mdj.Trajectory(grouped_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

# grouped_system.save_pdb("grouped_tryp_ben.pdb")

# finally we can combine these functions into something that recenters
# our box around a pair of groups (molecules for instance)

recentered_coords = recenter_pair(uncentered_system.xyz[0], box_lengths,
                                    prot_idxs, lig_idxs)

# write them out so we can visualize
recentered_system = mdj.Trajectory(recentered_coords, top,
                                unitcell_lengths=uncentered_system.unitcell_lengths,
                                unitcell_angles=uncentered_system.unitcell_angles)

recentered_system.save_pdb("recentered_tryp_ben.pdb")

