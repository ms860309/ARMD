import pymongo
from pymongo import MongoClient
import shutil
import os
import sys
from openbabel import pybel
from openbabel import openbabel as ob


class Connector(object):

    def __init__(self):
        #self.host = host
        #self.port = port
        self.server = 'mongodb+srv://jianyi:aa123@cluster0-wo5fn.gcp.mongodb.net/test?retryWrites=true&w=majority'
        #self.server = 'mongodb://localhost:27017/'
        #self.mongo_db = mongo_db
        self.client = self.connect()
        self.db = self.client['SnBEA_complete']

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=60000)
        return client


client = Connector()
db = getattr(Connector(), 'db')


def xyz_to_pyMol(xyz, cluster_bond_path=None):
    mol = next(pybel.readfile('xyz', xyz))
    if cluster_bond_path:
        m = pybel.ob.OBMol()
        m.BeginModify()
        for atom in mol:
            coords = [coord for coord in atom.coords]
            atomno = atom.atomicnum
            obatom = ob.OBAtom()
            obatom.thisown = 0
            obatom.SetAtomicNum(atomno)
            obatom.SetVector(*coords)
            m.AddAtom(obatom)
            del obatom

        with open(cluster_bond_path, 'r') as f:
            lines = f.read()
        cluster_bond = eval(lines)
        bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())
                 for bond in pybel.ob.OBMolBondIter(mol.OBMol)]
        bonds.extend(cluster_bond)
        for bond in bonds:
            m.AddBond(bond[0], bond[1], bond[2])
        pybelmol = pybel.Molecule(m)
        return pybelmol
    else:
        return mol

# debug
"""
reactions_collection = db['reactions']
acceptable_condition = ['forward equal to reactant and reverse equal to product',
                        'reverse equal to reactant and forward equal to product',
                        'forward equal to reactant but reverse does not equal to product',
                        'reverse equal to reactant but forward does not equal to product']

query = {'$and': 
                [
                {"unique":
                    {"$in":
                    ['new one']}},
                {'irc_equal':
                    {'$in':acceptable_condition}}
                ]
            }

reactions = list(reactions_collection.find(query))

for target in reactions:
    print('-----------')
    print('reactant_inchikey:{}'.format(target['reactant_inchi_key']))
    print('reactant smiles:{}'.format(target['reactant_smi']))
    print('product_inchikey:{}'.format(target['product_inchi_key']))
    print('product smiles:{}'.format(target['product_smi']))
    print('Barrier:{}'.format(target['barrier']))
    try:
        print('Delta H:{}'.format(target['delta_H']))
    except:
        print('Delta H:{}'.format(0))
    print('Generations:{}'.format(target['generations']))
    print('-----------')

"""
"""
qm_collection = db['qm_calculate_center']
statistics_collection = db['statistics']
#query = {'low_opt_status':"job_success"}
#a = list(qm_collection.find(query))
finished_reactant_list = []
for i in statistics_collection.find({}, {"_id": 0, "reactant_smiles": 0, "add how many products": 0, "generations": 0}):
    finished_reactant_list.append(i['reactant_inchi_key'])
    print(i['reactant_inchi_key'])
"""

"""
qm_collection = db['qm_calculate_center']
query = [{'$match':{'reactant_inchi_key':'OWCQMKVAAHGRRF-UHFFFAOYSA-N'}},
            {'$group':{'_id':'$reactant_inchi_key', 'reactant_mopac_hf':{'$min':'$reactant_mopac_hf'}}}]
a = list(qm_collection.aggregate(query))[0]['reactant_mopac_hf']
print(a)
"""

"""
qm_collection = db['qm_calculate_center']
query = {'ssm_status':'job_success'}
targets = list(qm_collection.find(query))

# use the checker.py path as the reference
checker_path = os.path.realpath(sys.argv[0])
ard_path = os.path.dirname(os.path.dirname(checker_path))
cluster_bond_path = os.path.join(ard_path, 'script/bonds.txt')

for target in targets:
    xyz = os.path.join(target['path'], 'ssm_product.xyz')
    pyMol_1 = xyz_to_pyMol(xyz, cluster_bond_path)
    prod_inchi_key = pyMol_1.write('inchiKey').strip()
    prod_smi = pyMol_1.write('can').split()[0]
    update_field = {'product_inchi_key':prod_inchi_key, 'product_smiles':prod_smi}
    qm_collection.update_one(target, {"$set": update_field}, True)

"""
"""
cluster_bond_path= '/home/jianyi/AutomaticReactionDiscovery/script/bonds.txt'
qm_collection = db['qm_calculate_center']
reactions_collection = db['reactions']
query = {"irc_opt_status":
                       {"$in":
                        ['job_success']}
        }
a = list(qm_collection.find({}))
b = list(reactions_collection.find({}))
count = 0
for i in a:
    try:
        reactant_path = os.path.join(i['path'], 'reactant.xyz')
        pyMol_3 = xyz_to_pyMol(reactant_path, cluster_bond_path=cluster_bond_path)
        reactant_inchi_key = pyMol_3.write('inchiKey').strip()
        reactant_smiles = pyMol_3.write('can').split()[0]
        qm_collection.update_one(i, {"$set": {'reactant_smiles':reactant_smiles}}, True)
    except:
        pass

for i in b:
    try:
        reactant_path = os.path.join(i['path'], 'reactant.xyz')
        pyMol_3 = xyz_to_pyMol(reactant_path, cluster_bond_path=cluster_bond_path)
        #reactant_inchi_key = pyMol_3.write('inchiKey').strip()
        reactant_smiles = pyMol_3.write('can').split()[0]
        qm_collection.update_one(i, {"$set": {'reactant_smiles':reactant_smiles}}, True)
    except:
        pass
"""
qm_collection = db['qm_calculate_center']
reactions_collection = db['reactions']
reactions = list(reactions_collection.aggregate([{
            '$group': {
                '_id': "$reaction",
                'barrier': {'$min': "$barrier"}
            }}
        ]))

for i in reactions:
    dir_path = list(reactions_collection.find({'reaction': i['_id'], 'barrier': i['barrier']}))[0]
    try:
        qmmm = dir_path['qmmm']
        if qmmm == 'Already insert':
            continue
    except:
        pass
    ard_qm_target = list(qm_collection.find({'path': dir_path['path']}))[0]
    # Prevent reactant equal to product but with different active site (maybe proton at the different oxygen)
    reactant_smiles = ard_qm_target['reactant_smiles'].split('.')
    product_smiles = ard_qm_target['product_smiles'].split('.')
    if len(reactant_smiles) > 1:
        reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs and 'C' in rs])
    else:
        reactant_part_smiles = set(reactant_smiles)

    if len(product_smiles) > 1:
        product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps and 'C' in ps])
    else:
        product_part_smiles = set(product_smiles)

    if reactant_part_smiles == product_part_smiles:
        continue

    qm_collection.update_one(ard_qm_target, {"$set": {"qmmm_opt_status": "job_unrun", "qmmm_freq_ts_status": "job_unrun"}}, True)
    reactions_collection.update_one(dir_path, {"$set": {'qmmm': 'Already insert'}}, True)