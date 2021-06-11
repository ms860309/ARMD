#!/usr/bin/env python
# -*- coding: utf-8 -*-


if __name__ == '__main__':

    # standard library imports
    import sys
    import os
    from os import path
    sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'code/ard'))
    sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'code/mol_graph'))
    sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'database'))

    # third party
    import argparse
    import logging
    from rdkit import Chem
    from rdkit import RDLogger
    from openbabel import pybel

    # local application imports
    from main import ARD, readInput, readXYZ, extract_bonds, extract_fixed_atoms_index

    # disable log
    RDLogger.DisableLog('rdApp.*')
    rootlogger = logging.getLogger()
    rootlogger.setLevel(logging.CRITICAL)
    pybel.ob.obErrorLog.SetOutputLevel(0)

    # Set up parser for reading the input filename from the command line
    parser = argparse.ArgumentParser(
        description='Automatic Reaction Mechanism Discovery')
    parser.add_argument('file', type=str, metavar='infile',
                        help='An input file describing the job options')
    parser.add_argument('reactant', default='reactant.xyz', type=str,
                        metavar='infile', help='An reactant xyz input file')
    parser.add_argument('-bonds', default='bonds.txt', type=str,
                        help='Manual specify bonds, (If you want to manually add bonds, set manual_bonds as 0 and manual_cluster_bond as 1 then put the bond you want to add in bonds.txt)', required=False)
    parser.add_argument('-fixed_atoms', default='fixed_atoms.txt', type=str,
                        help='Manual specify fixed atom index (start from 0)', required=False)
    parser.add_argument('-generations', default=1, type=int,
                        help='The network generation index', required=False)
    args = parser.parse_args()

    # Read input file
    input_file = path.abspath(args.file)
    ard_path = path.dirname(path.abspath(args.file))
    reactant_file = path.abspath(args.reactant)
    kwargs = readInput(input_file)

    if kwargs['fixed_atoms'] == '1':
        index = extract_fixed_atoms_index(args.fixed_atoms)
        kwargs['fixed_atoms'] = index
    else:
        kwargs['fixed_atoms'] = None

    # Manual set up bonds
    if kwargs['manual_bonds'] == '1':
        bonds = extract_bonds(args.bonds)
        kwargs['bonds'] = bonds
    else:
        if kwargs['manual_cluster_bond'] == '1':
            bonds = extract_bonds(args.bonds)
            kwargs['manual_cluster_bond'] = bonds
        else:
            kwargs['manual_cluster_bond'] = None
        kwargs['bonds'] = None
    OBMol, reactant_graph = readXYZ(reactant_file, kwargs['bonds'], kwargs['manual_cluster_bond'], kwargs['fixed_atoms'])

    kwargs['reactant'] = OBMol
    kwargs['graph'] = reactant_graph

    # Set output directory
    output_dir = path.abspath(path.dirname(input_file))
    kwargs['output_dir'] = output_dir
    kwargs['generations'] = args.generations
    kwargs['ard_path'] = ard_path
    kwargs['config_path'] = path.join(path.dirname(ard_path), 'config')
    kwargs['reactant_path'] = reactant_file

    # Execute job
    ard = ARD(**kwargs)
    ard.executeXYZ(**kwargs)
