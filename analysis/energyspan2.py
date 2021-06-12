import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def remove_same_element(lst):
    seen = set()
    identical = []
    for element in lst:
        if element not in seen:
            seen.add(element)
            identical.append(element)
    return identical


def update_energy():

    # forward reaction
    rxn_energy_dict.update({'{}_{}'.format(x, y): rxn_energy})
    barrier_dict.update({'{}_{}'.format(x, y): barrier})
    # backward reaction
    if '{}_{}'.format(x, y) not in rxn_edge:
        rxn_edge.append('{}_{}'.format(y, x))
    rxn_energy_dict.update({'{}_{}'.format(y, x): -rxn_energy})
    barrier_dict.update({'{}_{}'.format(y, x): barrier - rxn_energy})


def generate_energy_profile(node_list):
    profile = [0.0]
    for index, node in enumerate(node_list):
        if index > 0:
            edge = '{}_{}'.format(node_list[index - 1], node_list[index])
            if edge in rxn_energy_dict.keys():
                barrier = barrier_dict[edge] + profile[-1]
                state = rxn_energy_dict[edge] + profile[-1]
                profile.append(barrier)
                profile.append(state)
            else:
                return 'missing reaction'
                # return 'path: {} is not recommend: {} reaction not found.'.format(node_list, edge)
        if node == node_list[-1] and len(profile) == 2 * len(node_list) - 1:
            # return 'possible path:{}, profile:{}'.format(node_list, profile)
            return node_list, profile
            # return profile


def energy_span(energy_list):
    energy_profile = energy_list
    if energy_profile[0] != 0:
        print("Set initial state to reference state!")

    energy_ts = energy_profile[1::2]
    energy_int = [energy for energy in energy_profile if energy not in energy_ts]
    energy_rxn = energy_profile[-1] - energy_profile[0]
    energy_span_matrix = np.zeros((len(energy_int), len(energy_ts)))
    for j, ts in enumerate(energy_ts):
        for index, intermediate in enumerate(energy_int):
            if energy_profile.index(ts) < energy_profile.index(intermediate):
                energy_span_matrix[index, j] = ts - intermediate + energy_rxn
            else:
                energy_span_matrix[index, j] = ts - intermediate
    span = np.amax(energy_span_matrix)
    return span


def output_form(start, end, index, eff):
    repeat_edge = []
    # propose path with the (second) minimum energy span
    if start != '0' or end in usr_end:
        f.write('****user specified****\n')
    f.write('from {} to {}: \n'.format(start, end))
    f.write('path: {}\n'.format(', '.join(sub_path_node_list[index])))
    f.write('energy profile (kJ/mol): {}\n'.format(', '.join(map(str, sub_eng_profile[index]))))
    f.write('energy span (kJ/mol): {}\n'.format(eff))
    f.write('-------------\n')


current_dir = os.getcwd()
input_network_dir = os.path.join(current_dir, 'combine_test/network.txt')
output_label_dir = os.path.join(current_dir, 'combine_test/label.txt')
refine_dir = os.path.join(current_dir, 'combine_test/network-refine.txt')
output_graph_dir = os.path.join(current_dir, 'combine_test/graph.png')
output_energy_span_dir = os.path.join(current_dir, 'combine_test/output.txt')


with open(input_network_dir, "r") as original_file:
    rough_content = original_file.readlines()

# identify all species in the input file
species_list = []
for i, line in enumerate(rough_content):
    if line.__contains__('reactant_inchi_key:') or line.__contains__('product_inchi_key:'):
        species_list.append(line.split(':')[1])
node_identical = remove_same_element(species_list)

# create a dictionary {species_inchi_key: node_label}
label_dict = {}
with open(output_label_dir, "w") as label:
    for i, node in enumerate(node_identical):
        label.write('{} {}'.format(i, node))
        label_dict.update({node.strip(): str(i)})

with open(refine_dir, "w") as out:
    for line in rough_content:
        # read replace the string and write to output file
        if line.__contains__('reactant_inchi_key:') or line.__contains__('product_inchi_key:'):
            key = line.split(':')[1].strip()
            out.write(line.replace(key, label_dict[key]))
        elif line.__contains__('barrier:') or line.__contains__('delta_H:'):
            out.write(line)
        elif line.__contains__('-----------------') or line.__contains__('generations:'):
            out.write(line)

with open(refine_dir, "r") as network_file:
    content = network_file.readlines()

G = nx.Graph()
rxn_edge, all_input_edge = [], []
rxn_energy_dict, barrier_dict = {}, {}
for i, line in enumerate(content):
    if line.__contains__('generations:'):
        x = content[i - 2].split(':')[1].strip()
        y = content[i - 1].split(':')[1].strip()
        rxn_energy = float(content[i + 2].split(':')[1].strip())
        barrier = float(content[i + 1].split(':')[1].strip())
        # all_input_edge = [edge, barrier, rxn_energy]
        all_input_edge.append(['{}_{}'.format(x, y), barrier, rxn_energy])

        if '{}_{}'.format(x, y) not in rxn_edge and '{}_{}'.format(y, x) not in rxn_edge:
            G.add_edge(x, y)
            update_energy()

        elif '{}_{}'.format(x, y) in rxn_edge:
            if barrier < barrier_dict['{}_{}'.format(x, y)]:
                update_energy()


all_input_edge = np.array(sorted(all_input_edge))
unique, counts = np.unique(all_input_edge[:, 0], return_counts=True)
frequencies = np.asarray((unique, counts)).T.tolist()
repeat = [x for x, y in frequencies if y != '1']

# shells = [['0'], ['1', '2', '3', '4', '5'], ['6', '7', '8', '9', '10'],
#           ['11','12','13','14','15','16'],['17','18','19','20','21'],['22','23']]

nx.draw(G, pos=nx.spring_layout(G), with_labels=True)
plt.draw()
# plt.show()
plt.savefig(output_graph_dir)

# for user specified node (or inchi_key ex. str(label_dict['BAYVFPVHQJVOKP-XBOSCTMBSA-N']))
usr_start = ['0']
usr_end = []

# for end nodes
end_node_list = [n for n, d in G.degree() if d == 1]  # degree sequence

f = open(output_energy_span_dir, "w")
for start_node in usr_start:
    for end_node in end_node_list + usr_end:
        possible_energy_span, sub_path_node_list, sub_eng_profile = [], [], []
        for path in nx.all_simple_paths(G, source=start_node, target=end_node):
            path_node_list = [node for node in path]

            if generate_energy_profile(path_node_list) != 'missing reaction':
                path_node_list, eng_profile = generate_energy_profile(path_node_list)
                sub_path_node_list.append(path_node_list)
                sub_eng_profile.append(eng_profile)
                possible_energy_span.append(energy_span(eng_profile))

        eff_span = min(possible_energy_span)
        min_index = possible_energy_span.index(eff_span)
        output_form(start_node, end_node, min_index, eff_span)

        # consider path with the second minimum energy span,
        # if the energy difference is within 5 kJ/mol
        if (sorted(possible_energy_span)[1] - sorted(possible_energy_span)[0]) <= 5:
            eff_span_2 = sorted(possible_energy_span)[1]
            min_index_2 = possible_energy_span.index(sorted(possible_energy_span)[1])
            f.write('****second min path****\n')
            output_form(start_node, end_node, min_index_2, eff_span_2)

f.close()

# consider adsorption/desorption state of catalyst
