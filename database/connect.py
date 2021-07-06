from ntpath import join
import pymongo
from pymongo import MongoClient
import shutil
import os
import sys
from openbabel import pybel
from openbabel import openbabel as ob
from extract_time import find_job_output_file, get_job_run_time

class Connector(object):

    def __init__(self):
        #self.host = host
        #self.port = port
        self.server = 'mongodb+srv://jianyi:aa123@cluster0-wo5fn.gcp.mongodb.net/test?retryWrites=true&w=majority'
        #self.server = 'mongodb://localhost:27017/'
        #self.mongo_db = mongo_db
        self.client = self.connect()
        self.db = self.client['Final']

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
"""
qm_collection = db['qm_calculate_center']
reactions_collection = db['reactions']
reactions = list(reactions_collection.find({'qmmm':'already insert to qm'}))

for reaction in reactions:
    qmmm_reactant_path = os.path.join(reaction['path'], 'QMMM_REACTANT')
    qmmm_product_path = os.path.join(reaction['path'], 'QMMM_PRODUCT')
    qmmm_ts_path = os.path.join(reaction['path'], 'QMMM_TS')
    shutil.rmtree(qmmm_reactant_path)
    shutil.rmtree(qmmm_product_path)
    shutil.rmtree(qmmm_ts_path)
    # qm_target = list(qm_collection.find({'path':reaction['path']}))[0]
    # reactions_collection.update_one(reaction, {"$unset": {'qmmm_ts_freq_status':"", 'qmmm_freq_reactant_status':"", 'qmmm_freq_product_status':"", 'qmmm_ts_energy':"", 'qmmm_freq_reactant_energy':'', 'qmmm_freq_product_energy':""}}, True)
    # qm_collection.update_one(qm_target, {"$unset": {'qmmm_freq_ts_jobid':'', 'qmmm_freq_opt_status':'', 'qmmm_freq_product_jobid':'', 'qmmm_freq_product_status':'', 'qmmm_freq_reactant_jobid':'', 'qmmm_freq_reactant_status':'', 'qmmm_ts_freq_status':'', 'qmmm_ts_freq_jobid':'', 'qmmm_freq_opt_product_status':'', 'qmmm_freq_opt_reactant_status':'', 'qmmm_freq_opt_reactant_jobid':"", 'qmmm_freq_opt_product_jobid':''}, "$set": {'qmmm_freq_ts_status':'job_unrun', 'qmmm_freq_opt_product_restart_times':0, 'qmmm_freq_opt_reactant_restart_times':0, 'qmmm_freq_ts_restart_times':0, 'qmmm_opt_status':'job_unrun'}}, True)
"""
"""
reactions_collection = db['reactions']
reactions = list(reactions_collection.find({}))
for reaction in reactions:
    try:
        qmmm_barrier = reaction['qmmm_delta_H']
    except:
        continue
    print(f"{qmmm_barrier}   {reaction['delta_H']}")
"""

"""
reactions_collection = db['reactions']
reactions = list(reactions_collection.find({}))
for reaction in reactions:
    try:
        if reaction['qmmm_sp_reactant_status'] == 'job_success' and reaction['qmmm_sp_product_status'] == 'job_success' and reaction['qmmm_sp_ts_status'] == 'job_success' and reaction['qmmm_unique'] == 'unique':
            data_path = os.path.join(os.getcwd(), 'data')
            qmmm_reactant_dirpath = os.path.join(reaction['path'], 'QMMM_REACTANT')
            qmmm_product_dirpath = os.path.join(reaction['path'], 'QMMM_PRODUCT')
            qmmm_ts_dirpath = os.path.join(reaction['path'], 'QMMM_TS')
            dir_name = os.path.basename(reaction['path'])
            #a = os.path.join(data_path, dir_name)
            #os.mkdir(a)
            #shutil.copytree(qmmm_reactant_dirpath, os.path.join(a, 'QMMM_REACTANT'))
            #shutil.copytree(qmmm_product_dirpath, os.path.join(a, 'QMMM_PRODUCT'))
            #shutil.copytree(qmmm_ts_dirpath, os.path.join(a, 'QMMM_TS'))
            print(f"reactant_inchi_key:{reaction['reactant_inchi_key']}")
            print(f"product_inchi_key:{reaction['product_inchi_key']}")
            print(f"reactant_smiles:{reaction['reactant_smiles']}")
            print(f"product_smiles:{reaction['product_smiles']}")
            print(f"generations:{reaction['generations']}")
            print(f"qmmm_barrier:{reaction['qmmm_barrier']}")
            print(f"qmmm_delta_H:{reaction['qmmm_delta_H']}")
            print(f"dir_name:{dir_name}")
            print('------------------')
    except:
        continue
"""

"""
reactions_collection = db['reactions']
reactions = list(reactions_collection.find({}))
ssm_barriers, irc_barriers, irc_delta_Hs, qmmm_barriers, xTB_delta_Hs = [], [], [], [], []
for reaction in reactions:
    qm_collection = db['qm_calculate_center']
    target = list(qm_collection.find({'path':reaction['path']}))[0]
    try:
        ssm_barrier = target['ssm_barrier']
    except:
        ssm_barrier = 0
    try:
        irc_barrier = target['barrier']
    except:
        irc_barrier = 0
    try:
        irc_delta_H = reaction['delta_H']
    except:
        irc_delta_H = 0
    try:
        qmmm_barrier = target['qmmm_barrier']
    except:
        qmmm_barrier = 0
    try:
        xTB_delta_H = (target['product_xtb_hf'] - target['reactant_xtb_hf'])*627.5095
    except:
        xTB_delta_H = 0
    ssm_barriers.append(ssm_barrier)
    irc_barriers.append(irc_barrier)
    irc_delta_Hs.append(irc_delta_H)
    qmmm_barriers.append(qmmm_barrier)
    xTB_delta_Hs.append(xTB_delta_H)
    print(f'ssm_barrier:{ssm_barrier}')
    print(f'irc_barrier:{irc_barrier}')
    print(f'qmmm_barrier:{qmmm_barrier}')
    print(f"xTB_delta_H:{xTB_delta_H}")
    print(f'irc_delta_H:{irc_delta_H}')
    print('-----')

print('------------------')
print('ssm_barrier irc_barrier irc_delta_H qmmm_barrier xTB_delta_H')
for i,j,k,l,m in zip(ssm_barriers, irc_barriers, irc_delta_Hs, qmmm_barriers, xTB_delta_Hs):
    print(f'{i} {j} {k} {l} {m}')
"""

# qm_collection = db['qm_calculate_center']
# generations = [1, 2, 3, 4, 5]
# for gen in generations:
#     targets = list(qm_collection.find({'generations':gen}))
#     t = 0
#     for target in targets:
#         try:
#             a = int(target['irc_run_time'])*8
#         except:
#             a = 0
#         t += a
#     print(t)
"""
reactions_collection = db['reactions']
reactions = list(reactions_collection.find({}))
parse = []
count = 0
for reaction in reactions:
    try:
        if reaction['qmmm_sp_reactant_status'] == 'job_success' and reaction['qmmm_sp_product_status'] == 'job_success' and reaction['qmmm_sp_ts_status'] == 'job_success':
            _reaction = reaction['reaction']
            new_reaction = [_reaction[1], _reaction[0]]
            query = {'reaction':new_reaction}
            targets = list(reactions_collection.find(query))
            if reaction['path'] not in parse:
                for target in targets:
                    try:
                        if target['qmmm_sp_reactant_status'] == 'job_success' and target['qmmm_sp_product_status'] == 'job_success' and target['qmmm_sp_ts_status'] == 'job_success':
                            if target['path'] not in parse:
                                if target['qmmm_barrier'] - abs(reaction['qmmm_delta_H'] - reaction['qmmm_barrier']) > 5.0:
                                    print(target['qmmm_barrier'])
                                    print((reaction['qmmm_delta_H'] - reaction['qmmm_barrier']))
                                    print(target['path'])
                                else:
                                    reactions_collection.update_one(target, {"$set": {'qmmm_unique':"duplicate"}}, True)
                    except:
                        continue
            parse.append(reaction['path'])
    except:
        continue
"""

a="""QVGWVKLSHOQQBD-UHFFFAOYSA-N_2
XXUNXJUVVOLMDV-INDDKUIDSA-N_2
MRRZPFVYKCXOHV-UHFFFAOYSA-N_3
XXUNXJUVVOLMDV-INDDKUIDSA-N_1
OEDKBVBJDCKSOD-NBINMXDDSA-N_8
OEDKBVBJDCKSOD-WUPBIDMLSA-N_8
QVGWVKLSHOQQBD-UHFFFAOYSA-N_13
XXUNXJUVVOLMDV-UWKOMEOWSA-N_3
OYKSNHGFUKZQIC-UHFFFAOYSA-M_2
VQGFTDSQVVSIND-UHFFFAOYSA-N_3
ZRVNRYZCMKRAOZ-WUPBIDMLSA-N_4
OEDKBVBJDCKSOD-WUPBIDMLSA-N_13
SXYPHYKCDMDNMD-NBINMXDDSA-N_12
QUOGADVAICELDI-UHFFFAOYSA-N_14
SXYPHYKCDMDNMD-WUPBIDMLSA-N_11
SXYPHYKCDMDNMD-WUPBIDMLSA-N_14
GXDXQKACIWXXRY-NBINMXDDSA-N_13
QUOGADVAICELDI-UHFFFAOYSA-N_29
QVEKDYUISLHVHM-UHFFFAOYSA-N_24
HCHLEFUYHPKDDW-UWKOMEOWSA-N_1
ANUGDPPERUOYDH-UHFFFAOYSA-N_4
QUOGADVAICELDI-UHFFFAOYSA-N_46
SXYPHYKCDMDNMD-WUPBIDMLSA-N_25
QZVZJWXRGGCVGZ-WYGWCQTDSA-N_18
GXDXQKACIWXXRY-NBINMXDDSA-N_25
XXUNXJUVVOLMDV-PJUHLLTHSA-N_22
QVEKDYUISLHVHM-UHFFFAOYSA-N_35
BAYVFPVHQJVOKP-XBOSCTMBSA-N_19
QVEKDYUISLHVHM-UHFFFAOYSA-N_37
GXDXQKACIWXXRY-WUPBIDMLSA-N_20
BAYVFPVHQJVOKP-JGMWHDQESA-N_9
BAYVFPVHQJVOKP-UHFFFAOYSA-N_6
QWEZHXPVOWMARL-UHFFFAOYSA-N_4
SXYPHYKCDMDNMD-WUPBIDMLSA-N_32
OEDKBVBJDCKSOD-NBINMXDDSA-N_49
QVEKDYUISLHVHM-UHFFFAOYSA-N_46
SXYPHYKCDMDNMD-WUPBIDMLSA-N_35
ANUGDPPERUOYDH-UHFFFAOYSA-N_12
AWMDTVKXNQMJAY-UHFFFAOYSA-N_13
PGKPLWKYRUGQFB-NBINMXDDSA-N_23
OEDKBVBJDCKSOD-NBINMXDDSA-N_54
IVAQFJDWHKIGCA-UHFFFAOYSA-N_11
QJVZLJGSAOBMKY-UHFFFAOYSA-N_99
QVGWVKLSHOQQBD-UHFFFAOYSA-N_97
OEDKBVBJDCKSOD-WUPBIDMLSA-N_63
ANUGDPPERUOYDH-UHFFFAOYSA-N_48
YHHLECDPHPBLAQ-UHFFFAOYSA-N_37
QJVZLJGSAOBMKY-UHFFFAOYSA-N_115
IJRSHDCJFOGPCF-UHFFFAOYSA-N_16
QVEKDYUISLHVHM-UHFFFAOYSA-N_62
AZOYPMAKCCIFLG-UHFFFAOYSA-N_14
ANUGDPPERUOYDH-UHFFFAOYSA-N_83
ZLULKZQZVFEZGH-UHFFFAOYSA-N_41
FWXBZNPNRVQROJ-WYGWCQTDSA-N_10
ANUGDPPERUOYDH-UHFFFAOYSA-N_80
XXUNXJUVVOLMDV-UWKOMEOWSA-N_35
XXUNXJUVVOLMDV-PJUHLLTHSA-N_37"""
a = a.split('\n')

# for idx, i in enumerate(a):
#     dir_name = f'{i}'
#     _name = os.path.join('/home/jianyi/final/AutomatedReactionMechanismDiscovery/database/data/', dir_name)
#     os.mkdir(_name)
#     os.mkdir(os.path.join(_name, 'reactant'))
#     os.mkdir(os.path.join(_name, 'product'))
#     os.mkdir(os.path.join(_name, 'ts'))
#     reaction_path = f'/home/jianyi/final/AutomatedReactionMechanismDiscovery/reactions/{i}'
#     qmmm_reactant_dirpath = os.path.join(reaction_path, 'QMMM_REACTANT')
#     qmmm_product_dirpath =os.path.join(reaction_path, 'QMMM_PRODUCT')
#     qmmm_ts_dirpath = os.path.join(reaction_path, 'QMMM_TS')
#     shutil.copyfile(os.path.join(qmmm_reactant_dirpath, 'qmmm_freq_opt.out'), os.path.join(_name, 'reactant/qmmm_freq_opt.out'))
#     shutil.copyfile(os.path.join(qmmm_reactant_dirpath, 'qmmm_sp.out'), os.path.join(_name, 'reactant/qmmm_sp.out'))
#     shutil.copyfile(os.path.join(qmmm_reactant_dirpath, 'qmmm_final.xyz'), os.path.join(_name, 'reactant/qmmm_final.xyz'))
#     shutil.copyfile(os.path.join(qmmm_product_dirpath, 'qmmm_freq_opt.out'), os.path.join(_name, 'product/qmmm_freq_opt.out'))
#     shutil.copyfile(os.path.join(qmmm_product_dirpath, 'qmmm_sp.out'), os.path.join(_name, 'product/qmmm_sp.out'))
#     shutil.copyfile(os.path.join(qmmm_product_dirpath, 'qmmm_final.xyz'), os.path.join(_name, 'product/qmmm_final.xyz'))
#     shutil.copyfile(os.path.join(qmmm_ts_dirpath, 'qmmm_freq_ts.out'), os.path.join(_name, 'ts/qmmm_freq_ts.out'))
#     shutil.copyfile(os.path.join(qmmm_ts_dirpath, 'qmmm_sp.out'), os.path.join(_name, 'ts/qmmm_sp.out'))
#     shutil.copyfile(os.path.join(qmmm_ts_dirpath, 'qmmm_final.xyz'), os.path.join(_name, 'ts/qmmm_final.xyz'))
