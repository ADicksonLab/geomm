import numpy as np

def rmsd_one_frame(traj, ref, idx=[]):
    # calculates rmsd between frames of two trajectories: traj and ref
    if len(idx) == 0:
        idx = range(len(traj))
    # check if traj and ref have the same number of atoms
    if len(traj[idx]) != len(ref[idx]):
        raise Exception('RMSD calculation must use the same number of atoms in each array')
        
    return np.sqrt(np.sum(np.square(traj[idx, :] - ref[idx, :]),
                          axis=(0,1))/idx.shape[0])

