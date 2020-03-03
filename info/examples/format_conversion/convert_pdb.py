
# make sure to install dependencies in 'conversion.requirements.txt'
import os
import os.path as osp
import shutil

import sqlite3
import pandas as pd
import mdtraj as mdj
from wepy.util.mdtraj import mdtraj_to_json_topology
from wepy.util.json_top import json_top_chain_df, json_top_residue_df, json_top_atom_df

traj = mdj.load('../lysozyme_pxylene.pdb')

if osp.exists("outputs"):
    shutil.rmtree("outputs")
os.makedirs("outputs")


# the JSON format
json_top = mdtraj_to_json_topology(traj.top)

with open('outputs/lysozyme_pxylene.top.json', 'w') as wf:
    wf.write(json_top)

# FASTA residue sequence
fasta_str = traj.top.to_fasta(chain=0)
with open('outputs/lysozyme_pxylene.res.fasta', 'w') as wf:
    wf.write(fasta_str)

# topology tables

# you can get a table using mdtraj
mdj_atoms_df, bonds = traj.top.to_dataframe()

# but wepy provides the ability to get normalized versions for each
# level
chain_df = json_top_chain_df(json_top)
residue_df = json_top_residue_df(json_top)
atom_df = json_top_atom_df(json_top)

chain_df.to_csv('outputs/lysozyme_pxylene.chain.csv', index=False)
residue_df.to_csv('outputs/lysozyme_pxylene.residue.csv', index=False)
atom_df.to_csv('outputs/lysozyme_pxylene.atom.csv', index=False)


# to an SQLite3 database
db = sqlite3.Connection("outputs/lysozyme_pxylene.sqlite3")

chain_df.to_sql('chains', db)
residue_df.to_sql('residues', db)
atom_df.to_sql('atoms', db)

# to an HDF5 file
chain_df.to_hdf('outputs/lysozyme_pxylene.top.h5', 'chains')
residue_df.to_hdf('outputs/lysozyme_pxylene.top.h5', 'residues')
atom_df.to_hdf('outputs/lysozyme_pxylene.top.h5', 'atoms')

# to an excel spreadsheet
with pd.ExcelWriter('outputs/lysozyme_pxylene.top.xlsx', mode='r+') as writer:

    chain_df.to_excel(writer, sheet_name='chains')
    residue_df.to_excel(writer, sheet_name='residues')
    atom_df.to_excel(writer, sheet_name='atoms')
