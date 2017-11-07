import geomm.recentering
import geomm.rmsd
import mdtraj as mdj
import time
import numpy as np

pdbfile = 'complex.pdb'
dcdfile = 'cycle37.dcd'

traj = mdj.load_dcd(dcdfile,top=pdbfile)
n_frame = len(traj)
pdb = mdj.load_pdb(pdbfile)

lig_idxs = pdb.topology.select('resname "GST"')
binding_selection_idxs = [5,12,19,26,33,40,47,54]

box_lengths = np.array([np.array([4.07198,4.05594,4.04544]) for i in range(len(traj))])

start = time.time()
newpos = geomm.recentering.recenter_receptor_ligand(traj.xyz,box_lengths,ligand_idxs=lig_idxs,receptor_idxs=binding_selection_idxs)

traj_rec = mdj.Trajectory(newpos,pdb.top)
# superpose everything to frame zero
traj_rec.superpose(traj_rec,atom_indices=binding_selection_idxs)

d = np.zeros((n_frame,n_frame))
for i in range(n_frame-1):
    d[i][i] = 0
    for j in range(i+1,n_frame):
        d[i][j] = geomm.rmsd.rmsd_one_frame(traj_rec.xyz[i],traj_rec.xyz[j],lig_idxs)
        d[j][i] = d[i][j]

end = time.time()

np.savetxt("dmat.dat",d)
print("Time of calculation:",end-start)
