from openbabel import pybel
import networkx as nx
import matplotlib.pyplot as plt
import os
import os.path as path
import sys
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'database'))
from connect import db
import itertools
import math

def extract_data(barrier_threshold=200.0):

    reaction_collection = db['reactions']
    qm_collection = db['qm_calculate_center']
    reactions = list(reaction_collection.find({}))

    # Remove the reactant equal to product but with the different active site.
    # e.g. The Proton is at the different oxygen
    reactions_copy = reactions[:]
    for reaction in reactions_copy:
        same = False
        reactant_smiles = reaction['reactant_smiles'].split('.')
        product_smiles = reaction['product_smiles'].split('.')
        if len(reactant_smiles) > 1:
            reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs and 'C' in rs])
        else:
            reactant_part_smiles = set(reactant_smiles)

        if len(product_smiles) > 1:
            product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps and 'C' in ps])
        else:
            product_part_smiles = set(product_smiles)

        if reactant_part_smiles == product_part_smiles:
            same = True

        if same:
            reactions.remove(reaction)

        if reaction['manual_check'] != 'need check':
            reactions.remove(reaction)

    reactant_smi, product_smi, barrier, generations, er_smi, ep_smi = [], [], [], [], [], []
    for target in reactions:
        reactant_smiles = target['reactant_smiles'].split('.')
        product_smiles = target['product_smiles'].split('.')
        if len(reactant_smiles) > 1:
            reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs])
            reactant_part_smiles = '.'.join(reactant_part_smiles)
        else:
            reactant_part_smiles = list(set(reactant_smiles))[0]

        if len(product_smiles) > 1:
            product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps])
            product_part_smiles = '.'.join(product_part_smiles)
        else:
            product_part_smiles = list(set(product_smiles))[0]

        reactant_smi.append(reactant_part_smiles)
        product_smi.append(product_part_smiles)
        er_smi.append(target['reactant_smiles'])
        ep_smi.append(target['product_smiles'])
        # print(target['barrier'])
        barrier.append(round(target['barrier'], 2))
        generations.append(target['generations'])

        test_target = list(qm_collection.find({'path': target['path']}))[0]
        # dH = (test_target['product_xtb_hf'] - test_target['reactant_xtb_hf']) * 627.5095
    zipped = zip(reactant_smi, product_smi, generations, barrier, er_smi, ep_smi)

    return zipped

def draw(target_product = None):
    G = nx.DiGraph()  # create object
    eG = nx.DiGraph() # For get energy
    zipped = extract_data()
    _dict = {}
    labels = {}
    for i, j, l, k, m, n in list(zipped):
        if l == 1:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='r')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        elif l == 2:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='g')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        elif l == 3:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='b')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        elif l == 4:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='y')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        elif l == 5:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='c')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        elif l == 6:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='m')
            eG.add_edge(m, n)
            _dict[(i, j)] = k
        else:
            if G.has_edge(i, j):
                continue
            G.add_edge(i, j, color='k')
            eG.add_edge(m, n)
            _dict[(i, j)] = k

    plt.figure(figsize=(10, 10))
    
    colors = nx.get_edge_attributes(G, 'color').values()
    weights = nx.get_edge_attributes(G, 'weight').values()

    # pos = nx.circular_layout(G)
    # pos = nx.shell_layout(G)
    # pos = nx.spring_layout(G)
    pos = nx.kamada_kawai_layout(G, scale=3)

    nx.draw(G, pos,
            edge_color=colors,
            with_labels=True,
            node_color='green',
            font_size=5, 
            node_size=1500,
            connectionstyle='arc3, rad = 0.15')

    nx.draw_networkx_edge_labels(G, pos, edge_labels=_dict, font_size=8)
    nx.draw_networkx_labels(G,pos,labels,font_size=10,font_color='r')
    root_to_leaf_paths(eG, target_product = target_product)
    
    plt.savefig("simple_path_3.png", dpi=1000) # save as png
    plt.show()


def root_to_leaf_paths(G, target_product = None):
    reaction_collection = db['reactions']
    _nodes = []
    for node in G.nodes:
        paths = []
        energy_profiles = []
        if G.degree[node] == 1 and node not in _nodes:
            if target_product:
                for path in nx.all_simple_paths(G, source='CC(=O)CO', target=target_product):
                    paths.append(path)
                    _nodes.append(node)
            else:
                for path in nx.all_simple_paths(G, source='CC(=O)CO.[SiH3]O[Sn]1(O[SiH3])(O[SiH3])[OH][SiH3][OH]1', target=node):
                    paths.append(path)
                    _nodes.append(node)
            paths.sort()
            paths = list(k for k,_ in itertools.groupby(paths))
            for path in paths:
                energy_profile = [0]
                for i in range(len(path)):
                    try:
                        reactant_smiles = path[i]
                        product_smiles = path[i+1]
                        reactions = list(reaction_collection.aggregate([{
                                                                '$match':{
                                                                    'reactant_smiles':reactant_smiles,
                                                                    'product_smiles':product_smiles
                                                                }},{
                                                                '$group':{
                                                                        '_id': "$reaction",
                                                                        'barrier': {'$min': "$barrier"}
                                                                        }}
                                                                ]))[0]
                    except:
                        continue
                    target = list(reaction_collection.find({'reaction':reactions['_id'], 'barrier':reactions['barrier']}))[0]
                    if i == 0:
                        energy_profile.append(round(target['barrier'], 2))
                        energy_profile.append(round(target['delta_H'], 2))
                    else:
                        energy_profile.append(round(energy_profile[-1]+target['barrier'], 2))
                        energy_profile.append(round(energy_profile[-2]+target['delta_H'], 2))
                energy_profiles.append(energy_profile)
            if energy_profiles == []:
                continue
            print_energy_profile_inf(paths, energy_profiles)


def print_energy_profile_inf(paths, energy_profiles):
    min_barrier_energy = max(energy_profiles[0])
    min_idx = 0
    for idx, energy_profile in enumerate(energy_profiles):
        if max(energy_profile) < min_barrier_energy:
            min_barrier_energy = max(energy_profile)
            min_idx = idx
    
    significant_step = list(filter(
        lambda step: all([step in path for path in paths]), paths[list(map(len, paths)).index((min(map(len, paths))))]
        ))

    new_path = []
    for old_path in paths[min_idx]:
        tmp = []
        old_paths = old_path.split('.')
        for opath in old_paths:
            if 'Sn' not in opath and len(old_paths) > 1:
                tmp.append(opath)
            elif len(old_paths) == 1:
                tmp.append(opath)
        smiles = '.'.join(tmp)
        new_path.append(smiles)

    print('--------------------')
    print('Reaction profile:\n{}'.format(new_path))
    print('Energy profile:\n{}'.format(energy_profiles[min_idx]))
    print('Maximum barrier:\n{}'.format(min_barrier_energy))
    print('The significant steps of the reaction:\n{}'.format(significant_step))
    print('--------------------\n\n')


draw(target_product = None)
