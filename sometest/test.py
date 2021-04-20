import pybel
import openbabel as ob
import numpy as np
import os
import pymongo
from pymongo import MongoClient

def check_bond_length(coords, add_bonds):
    """
    Use reactant coordinate to check if the add bonds's bond length is too long.
    Return a 'list of distance'.
    """
    dist = []
    for bond in add_bonds:
        coord_vect_1 = coords[0][bond[0]]
        coord_vect_2 = coords[0][bond[1]]
        diff = coord_vect_1 - coord_vect_2
        dist.append(np.linalg.norm(diff))
    return dist

server = 'mongodb+srv://jianyi:aa123@cluster0-wo5fn.gcp.mongodb.net/test?retryWrites=true&w=majority'
client = MongoClient(server, serverSelectionTimeoutMS=2000)
db = client['network']
qm_collection = db['qm_calculate_center']

xyz = '/mnt/d/reactant.xyz'
OBMol = next(pybel.readfile('xyz', xyz))

atoms = tuple(atom.atomicnum for atom in OBMol)
coords = [atom.coords for atom in OBMol]
reactant_coords = [np.array(coords).reshape(len(atoms), 3)]

reaction_dir = '/mnt/d/reactions'
remote_path = '/home/jianyi/AutomaticReactionDiscovery/reactions/'
a = os.listdir(reaction_dir)
for i in a:
    real_path = os.path.join(reaction_dir, i)
    remote = os.path.join(remote_path, i)
    query = {'path':remote}
    targets = list(qm_collection.find(query))
    isomer_path = os.path.join(real_path, 'add_bonds.txt')
    with open(isomer_path, 'r') as f:
        lines = f.read().splitlines()
    tmp = []
    for line in lines:
        line = line.split()
        if line[0] == 'ADD':
            tmp.append((int(line[1])-1, int(line[2])-1))
    dist = check_bond_length(reactant_coords, tmp)
    with open('debug.txt', 'a') as f:
        f.write(i)
        f.write('\n')
        f.write(str(dist))
        f.write('\n')
        f.write('ssm_status:{}'.format(targets[0]['ssm_status']))
        f.write('\n')
        try:
            f.write('ts_status:{}'.format(targets[0]['ts_status']))
        except:
            f.write('ts_status:{}'.format('job_unrun'))
        f.write('\n')
        f.write('---------')
        f.write('\n')
