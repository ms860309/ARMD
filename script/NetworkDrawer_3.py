from openbabel import pybel
import networkx as nx
import matplotlib.pyplot as plt
import os
import os.path as path
import sys
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'database'))
from connect import db
import itertools

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

    reactant_smi, product_smi, barrier, generations = [], [], [], []
    for target in reactions:
        reactant_smiles = target['reactant_smiles'].split('.')
        product_smiles = target['product_smiles'].split('.')
        if len(reactant_smiles) > 1:
            reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs])
            reactant_part_smiles = '.'.join(reactant_part_smiles)
        else:
            reactant_part_smiles = set(reactant_smiles)[0]

        if len(product_smiles) > 1:
            product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps])
            product_part_smiles = '.'.join(product_part_smiles)
        else:
            product_part_smiles = set(product_smiles)[0]

        reactant_smi.append(reactant_part_smiles)
        product_smi.append(product_part_smiles)
        # print(target['barrier'])
        barrier.append(round(target['barrier'], 2))
        generations.append(target['generations'])

        test_target = list(qm_collection.find({'path': target['path']}))[0]
        # dH = (test_target['product_xtb_hf'] - test_target['reactant_xtb_hf']) * 627.5095
    zipped = zip(reactant_smi, product_smi, generations, barrier)

    return zipped

def draw():
    G = nx.DiGraph()  # create object

    zipped = extract_data()
    _dict = {}
    labels = {}
    for i, j, l, k in list(zipped):
        if l == 1:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='r')
            # _dict[(i, j)] = k
        elif l == 2:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='g')
            # _dict[(i, j)] = k
        elif l == 3:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='b')
            # _dict[(i, j)] = k
        elif l == 4:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='y')
            # _dict[(i, j)] = k
        elif l == 5:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='c')
            # _dict[(i, j)] = k
        elif l == 6:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='m')
            # _dict[(i, j)] = k
        else:
            if G.has_edge(i, j) :
                continue
            G.add_edge(i, j, color='k')
            # _dict[(i, j)] = k

    plt.figure(figsize=(8, 8))
    
    colors = nx.get_edge_attributes(G, 'color').values()
    weights = nx.get_edge_attributes(G, 'weight').values()

    #pos = nx.circular_layout(G)
    #pos = nx.shell_layout(G)
    pos = nx.spring_layout(G)
    # pos = nx.kamada_kawai_layout(G)
    nx.draw(G, pos,
            edge_color=colors,
            with_labels=True,
            node_color='white',
            font_size=8)
    
    #nx.draw_networkx_edge_labels(G, pos, font_size=8)
    nx.draw_networkx_labels(G,pos,labels,font_size=10,font_color='r')
    #root_to_leaf_paths(G)

    plt.savefig("simple_path_2.png")  # save as png
    plt.show()



draw()
