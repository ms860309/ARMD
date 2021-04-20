# standard library imports
import sys
import os
from os import path
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'code/mol_graph'))

# local application imports
from input_output import xyz_file_to_atoms
from species import Species
from graph import make_graph
from bond_rearrangement import get_bond_rearrangs



reactant = Species(xyz_file_to_atoms('../script/reactant.xyz'))
product = Species(xyz_file_to_atoms('../script/reactant.xyz'))

bonds = '../script/bonds.txt'
with open(bonds, 'r') as f:
    lines = f.read()
reactant_bonds = [(i[0]-1, i[1]-1) for i in eval(lines)]
product_bonds = [(1,2,1),(1,3,1),(1,4,1),(1,5,1),(1,13,1),(2,3,1),(2,4,1),(2,6,1),(3,5,1),(3,6,1),(3,10,1),(4,5,1),(4,6,1),(5,6,1),(7,14,1),(8,15,1),(9,10,1),(9,13,1),(10,11,1),(10,16,1),(11,12,1),(11,17,1),(11,20,1),(12,13,1),(12,18,1),(12,21,1),(13,14,1),(14,15,1),(14,19,1)]
product_bonds = [(i[0]-1, i[1]-1) for i in product_bonds]

make_graph(reactant, bond_list= reactant_bonds)
make_graph(product, bond_list= product_bonds)

a = get_bond_rearrangs(reactant, product)