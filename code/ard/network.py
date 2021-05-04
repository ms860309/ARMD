# standard library imports
import os
import time
import psutil
import shutil

# third party
from openbabel import pybel
from openbabel import openbabel
try:
    from rmgpy import settings
    from rmgpy.data.thermo import ThermoDatabase
    from rmgpy.molecule import Molecule
except:
    pass
from subprocess import Popen, PIPE
import numpy as np

# local application imports
import gen3D
import util
from quantum import QuantumError
from node import Node
from pgen import Generate
from mopac import Mopac
import mopac
from _xtb import XTB

# database
from connect import db

info = psutil.virtual_memory()


class Network(object):
    def __init__(self, logger, **kwargs):
        self.logger = logger
        self.cluster_bond = kwargs['manual_cluster_bond']
        self.forcefield = kwargs['forcefield']
        self.constraintff_alg = kwargs['constraintff_alg']
        self.mopac_method = kwargs['mopac_method']
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.ard_path = kwargs['ard_path']
        self.generations = kwargs['generations']
        self.method = kwargs["dh_cutoff_method"]
        self.constraint = kwargs['constraint_index']
        self.count = 0
        self.fixed_atom = kwargs['fixed_atom']
        self.use_irc = kwargs['use_irc']
        self.reactant_path = os.path.dirname(kwargs['reactant_path'])

    def genNetwork(self, mol_object, **kwargs):
        """
        Execute the automatic reaction discovery procedure.
        """
        # Database
        qm_collection = db['qm_calculate_center']
        config_collection = db['config']
        statistics_collection = db['statistics']
        targets = list(config_collection.find({'generations': 1}))
        config_collection.update_one(targets[0], {"$set": {'config_path':kwargs['config_path']}}, True)
        # Reactant information
        reactant_inchi_key = mol_object.write('inchiKey').strip()  # inchikey

        # Generate all possible products
        gen = Generate(mol_object, **kwargs)
        self.logger.info('Generating all possible products...')
        gen.generateProducts()
        prod_mols = gen.get_prods()
        add_bonds = gen.get_add_bonds()
        break_bonds = gen.get_break_bonds()
        prod_mols_filtered = []
        self.logger.info(f'{len(prod_mols)} possible products are generated\n')

        # Filter reactions based on standard heat of reaction  delta H
        if self.method.lower() == 'mopac':
            self.logger.info(f'Now use {self.method} to filter the delta H of reactions....\n')
            if self.generations == 1:
                os.mkdir(os.path.join(os.path.dirname(self.ard_path), 'reactions'))
                H298_reac = self.get_mopac_H298(mol_object)
                config_collection.update_one(targets[0], {"$set": {'reactant_energy': H298_reac, 'use_irc': self.use_irc}}, True)
                mol_object_copy = mol_object.copy()
                for prod_mol in prod_mols:
                    if self.filter_dh_mopac(mol_object, self.cluster_bond, prod_mol, add_bonds[prod_mols.index(prod_mol)], 
                                            break_bonds[prod_mols.index(prod_mol)], len(prod_mols), qm_collection, refH=None):
                        prod_mols_filtered.append(prod_mol)
                    # Recovery
                    mol_object_copy = mol_object.copy()
            else:
                H298_reac = targets[0]['reactant_energy']
                mol_object_copy = mol_object.copy()
                for prod_mol in prod_mols:
                    if self.filter_dh_mopac(mol_object, self.cluster_bond, prod_mol, add_bonds[prod_mols.index(prod_mol)], 
                                            break_bonds[prod_mols.index(prod_mol)], len(prod_mols), qm_collection, refH=H298_reac):
                        prod_mols_filtered.append(prod_mol)
                    # Recovery
                    mol_object_copy = mol_object.copy()
        elif self.method.lower() == 'xtb':
            self.logger.info('Now use {} to filter the delta H of reactions....\n'.format(self.method))
            if self.generations == 1:
                os.mkdir(os.path.join(os.path.dirname(self.ard_path), 'reactions'))
                H298_reac = self.get_xtb_H298(config_path=kwargs['config_path'])
                config_collection.update_one(targets[0], {"$set": {'reactant_energy': H298_reac, 'use_irc': self.use_irc}}, True)
                mol_object_copy = mol_object.copy()
                for prod_mol in prod_mols:
                    if self.filter_dh_xtb(mol_object, self.cluster_bond, prod_mol, add_bonds[prod_mols.index(prod_mol)], 
                                            break_bonds[prod_mols.index(prod_mol)], len(prod_mols), qm_collection, config_path = kwargs['config_path'], refH=None):
                        prod_mols_filtered.append(prod_mol)
                    mol_object.setCoordsFromMol(mol_object_copy)
            else:
                H298_reac = targets[0]['reactant_energy']
                mol_object_copy = mol_object.copy()
                for prod_mol in prod_mols:
                    if self.filter_dh_xtb(mol_object, self.cluster_bond, prod_mol, add_bonds[prod_mols.index(prod_mol)], 
                                            break_bonds[prod_mols.index(prod_mol)], len(prod_mols), qm_collection, config_path = kwargs['config_path'], refH=H298_reac):
                        prod_mols_filtered.append(prod_mol)
                    mol_object.setCoordsFromMol(mol_object_copy)
        else:
            self.logger.info('Now use {} to filter the delta H of reactions....\n'.format(self.method))
            # Load thermo database and choose which libraries to search
            thermo_db = ThermoDatabase()
            thermo_db.load(os.path.join(
                settings['database.directory'], 'thermo'))
            thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                      'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]
            if self.generations == 1:
                H298_reac = mol_object.getH298(thermo_db)
                update_field = {'reactant_energy': H298_reac}
                config_collection.update_one(targets[0], {"$set": update_field}, True)
            else:
                H298_reac = targets[0]['reactant_energy']
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_rmg(H298_reac, mol, thermo_db)]
            self.logger.info('Generate geometry........\n')
            for mol in prod_mols_filtered:
                index = prod_mols.index(mol)
                # Generate geometry and return path
                mol_object.gen3D(self.constraint, forcefield=self.forcefield, method=self.constraintff_alg, make3D=False)
                reac_mol_copy = mol_object.copy()
                dir_path = self.gen_geometry(mol_object, mol, reac_mol_copy, add_bonds[index], break_bonds[index])
                product_inchi_key = mol.write('inchiKey').strip()
                self.logger.info('\nReactant inchi key: {}\nProduct inchi key: {}\nReactant smiles: {}\nProduct smiles: {}\nDirectory path: {}\n'.format(
                                    reactant_inchi_key, product_inchi_key, mol_object.write('can').split(), mol.write('can').split(), dir_path))

                qm_collection.insert_one({
                    'reaction': [reactant_inchi_key, product_inchi_key],
                    'reactant_smiles': mol_object.write('can').split()[0],
                    'reactant_inchi_key': reactant_inchi_key,
                    'product_inchi_key': product_inchi_key,
                    'product_smiles': mol.write('can').split()[0],
                    'path': dir_path,
                    'ssm_status': 'job_unrun',
                    'generations': self.generations,
                    'use_irc': self.use_irc
                })

        self.logger.info('After delta H filter {} product remain.\n'.format(len(prod_mols_filtered)))
        # Generate geometry and insert to database
        statistics_collection.insert_one({
            'reactant_smiles': mol_object.write('can').split()[0],
            'reactant_inchi_key': reactant_inchi_key,
            'add how many products': len(prod_mols_filtered),
            'generations': self.generations})

    def filter_dh_rmg(self, H298_reac, prod_mol, thermo_db):
        """
        Filter threshold based on standard enthalpies of formation of reactants
        and products. Returns `True` if the heat of reaction is less than
        `self.dh_cutoff`, `False` otherwise.
        """
        H298_prod = prod_mol.getH298(thermo_db)
        dH = H298_prod - H298_reac
        if dH < self.dh_cutoff:
            return 1
        return 0

    def filter_dh_mopac(self, reac_obj, cluster_bond, prod_mol, form_bonds, break_bonds, total_prod_num, qm_collection, refH=None):
        self.count += 1
        mopac_object = Mopac(reac_obj, prod_mol, self.mopac_method, self.forcefield, self.constraintff_alg,
                             form_bonds, break_bonds, self.logger, total_prod_num, self.count, self.constraint, self.fixed_atom, cluster_bond)
        H298_reac, H298_prod = mopac_object.mopac_get_H298(self.reactant_path)

        if H298_prod == False or H298_reac == False:
            return 0
        self.logger.info(
            'Product energy calculate by mopac is {} kcal/mol and reactant is {} kcal/mol'.format(H298_prod, H298_reac))
        if refH:
            self.logger.info('In the {} generations, reactant hf use {} instead.'.format(self.generations, refH))
            dH = H298_prod - refH
        else:
            dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            self.logger.info('Delta H is {}, smaller than threshold'.format(dH))

            reactant_output = os.path.join(self.reactant_path, 'tmp/reactant.out')
            product_output = os.path.join(self.reactant_path, 'tmp/product.out')

            dir_path = self.mopac_output(reactant_output, product_output, form_bonds, break_bonds, prod_mol)
            reactant_inchi_key = reac_obj.write('inchiKey').strip()
            product_inchi_key = prod_mol.write('inchiKey').strip()
            self.logger.info('\nReactant inchi key: {}\nProduct inchi key: {}\nReactant smiles: {}\nProduct smiles: {}\nDirectory path: {}\n'.format(
                            reactant_inchi_key, product_inchi_key, reac_obj.write('can').split(), prod_mol.write('can').split(), dir_path))
            qm_collection.insert_one({
                'reaction': [reactant_inchi_key, product_inchi_key],
                'reactant_smiles': reac_obj.write('can').split()[0],
                'reactant_inchi_key': reactant_inchi_key,
                'product_inchi_key': product_inchi_key,
                'product_smiles': prod_mol.write('can').split()[0],
                'reactant_mopac_hf': H298_reac,
                'product_mopac_hf': H298_prod,
                'path': dir_path,
                'ssm_status': 'job_unrun',
                'generations': self.generations,
                'use_irc': self.use_irc
            })
            self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
            return 1
        else:
            self.logger.info('Delta H is {}, greater than threshold'.format(dH))
            self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
            return 0

    def filter_dh_xtb(self, reac_obj, cluster_bond, prod_mol, form_bonds, break_bonds, total_prod_num, qm_collection, config_path, refH=None):
        self.count += 1
        xtb_object = XTB(reac_obj, prod_mol, self.forcefield, self.constraintff_alg, form_bonds, break_bonds,
                         self.logger, total_prod_num, self.count, self.constraint, self.fixed_atom, cluster_bond)
        H298_reac, H298_prod = xtb_object.xtb_get_H298(self.reactant_path, config_path)

        if H298_prod == False or H298_reac == False:
            return 0
        self.logger.info('Product energy calculate by xtb is {} Eh and reactant is {} Eh'.format(H298_prod, H298_reac))
        if refH:
            self.logger.info('In the {} generations, reactant hf use {} instead.'.format(self.generations, refH))
            dH = (H298_prod - refH) * 627.5095  # convert Eh to kcal/mol
        else:
            dH = (H298_prod - H298_reac) * 627.5095

        if dH < self.dh_cutoff:
            self.logger.info('Delta H is {}, smaller than threshold'.format(dH))
            reactant_output = os.path.join(self.reactant_path, 'tmp/reactant.xyz')
            product_output = os.path.join(self.reactant_path, 'tmp/product.xyz')

            dir_path = self.xtb_output(reactant_output, product_output, form_bonds, break_bonds, prod_mol)
            reactant_inchi_key = reac_obj.write('inchiKey').strip()
            product_inchi_key = prod_mol.write('inchiKey').strip()
            self.logger.info('\nReactant inchi key: {}\nProduct inchi key: {}\nReactant smiles: {}\nProduct smiles: {}\nDirectory path: {}\n'.format(
                            reactant_inchi_key, product_inchi_key, reac_obj.write('can').split(), prod_mol.write('can').split(), dir_path))

            qm_collection.insert_one({
                'reaction': [reactant_inchi_key, product_inchi_key],
                'reactant_smiles': reac_obj.write('can').split()[0],
                'reactant_inchi_key': reactant_inchi_key,
                'product_inchi_key': product_inchi_key,
                'product_smiles': prod_mol.write('can').split()[0],
                'reactant_xtb_hf': H298_reac,
                'product_xtb_hf': H298_prod,
                'path': dir_path,
                'ssm_status': 'job_unrun',
                'generations': self.generations,
                'use_irc': self.use_irc
            })
            self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
            return 1
        else:
            self.logger.info('Delta H is {}, greater than threshold'.format(dH))
            self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
            return 0

    def get_mopac_H298(self, mol_object, charge=0, multiplicity='SINGLET'):
        tmpdir = os.path.join(self.reactant_path, 'tmp')
        reactant_path = os.path.join(tmpdir, 'reactant.mop')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)

        reac_geo = str(mol_object.toNode()).splitlines()
        reactant_geometry = []
        for idx, i in enumerate(reac_geo):
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            if not self.constraint:
                l = " 1 ".join(k)
            else:
                if idx in self.constraint:
                    l = " 0 ".join(k)
                else:
                    l = " 1 ".join(k)
            out = atom + l
            reactant_geometry.append(out)
        reactant_geometry = "\n".join(reactant_geometry)

        with open(reactant_path, 'w') as f:
            f.write("NOSYM CHARGE={} {} {}\n\n".format(
                charge, multiplicity, self.mopac_method))
            f.write("\n{}".format(reactant_geometry))
        mopac.runMopac(tmpdir, 'reactant.mop')
        mol_hf = mopac.getHeatofFormation(tmpdir, 'reactant.out')
        return float(mol_hf)

    def get_xtb_H298(self, config_path):
        tmpdir = os.path.join(self.reactant_path, 'tmp')
        reactant_path = os.path.join(tmpdir, 'reactant.xyz')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)

        shutil.copyfile(os.path.join(self.ard_path, 'reactant.xyz'), reactant_path)
        try:
            mol_hf = self.runXTB(tmpdir, config_path, 'reactant.xyz')
        except:
            raise Exception('The initial reactant energy calculation by xtb is fail.')
        return float(mol_hf)

    def unique_key_filterIsomorphic_itself(self, base):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        base_unique = [mol.write('inchiKey').strip() for mol in base]
        #base_unique = [mol.toRMGMolecule().to_inchi_key() for mol in base]
        result = [base[base_unique.index(i)] for i in set(base_unique)]
        return result

    def gen_geometry(self, reactant_mol, product_mol, reactant_mol_copy, add_bonds, break_bonds, **kwargs):
        # Initial optimization
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
        reactant_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        product_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        # Arrange
        # If arrange error can use try
        arrange3D = gen3D.Arrange3D(
            reactant_mol, product_mol, self.constraint, self.fixed_atom, self.cluster_bond)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            self.logger.info(msg)

        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        ff.Setup(Hatom.OBMol)
        reactant_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        ff.Setup(Hatom.OBMol)
        product_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        ff.Setup(Hatom.OBMol)

        reactant = reactant_mol.toNode()
        product = product_mol.toNode()
        self.logger.info(
            'Reactant and product geometry is :\n{}\n****\n{}'.format(str(reactant), str(product)))
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        b_dirname = product_mol.write('inchiKey').strip()
        dirname = self.dir_check(subdir, b_dirname)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir
        #self.makeInputFile(reactant, product, **kwargs)
        self.makeCalFile(reactant, 'reactant.xyz', **kwargs)
        self.makeCalFile(product, 'product.xyz', **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)

        reactant_mol.setCoordsFromMol(reactant_mol_copy)
        return output_dir

    def mopac_output(self, reactant_output, product_output, add_bonds, break_bonds, prod_mol, **kwargs):
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        b_dirname = prod_mol.write('inchiKey').strip()
        dirname = self.dir_check(subdir, b_dirname)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir

        rsymbol, rgeometry = self.get_mopac_opt_geometry(reactant_output)
        psymbol, pgeometry = self.get_mopac_opt_geometry(product_output)

        #self.makeInputFile(reactant, product, **kwargs)
        self.makeCalFile(reactant_output, 'reactant.xyz', symbol=rsymbol, geometry=rgeometry, **kwargs)
        self.makeCalFile(product_output, 'product.xyz', symbol=psymbol, geometry=pgeometry, **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)
        return output_dir

    def xtb_output(self, reactant_output, product_output, add_bonds, break_bonds, prod_mol, **kwargs):
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        
        b_dirname = prod_mol.write('inchiKey').strip()
        dirname = self.dir_check(subdir, b_dirname)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir

        shutil.copyfile(reactant_output, os.path.join(output_dir, 'reactant.xyz'))
        shutil.copyfile(product_output, os.path.join(output_dir, 'product.xyz'))
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)
        return output_dir

    @staticmethod
    def prepare_mol(reactant_mol):
        """
        Covert xyz to mol.
        """
        atom_symbol = []
        coords = []
        for mol in reactant_mol.mols:
            atoms = tuple(atom.atomicnum for atom in mol)
            for atom in atoms:
                atom_symbol.append(openbabel.GetSymbol(atom))
            coords.append(mol.toNode().coords)
        coords = np.concatenate(coords)
        coords = coords.reshape(-1, 3)
        coords = "\n".join(
                ["{} {} {} {}".format(a, *c) for a, c in zip(atom_symbol, coords)]
        )
        xyz = f"{len(atom_symbol)}\n\n{coords}"
        mol = pybel.readstring("xyz", xyz)
        return gen3D.Molecule(mol.OBMol)

    @staticmethod
    def dir_check(subdir, b_dirname, num = 1):
        """
        When parallely run job, the dir is constructed but data is not on database yet
        """
        check = True
        while check:
            new_name = '{}_{}'.format(b_dirname, num)
            if os.path.exists(os.path.join(subdir, new_name)):
                num += 1
            else:
                check = False
        return new_name

    @staticmethod
    def makeInputFile(reactant, product, **kwargs):
        """
        Create input file for TS search and return path to file.
        """
        output = os.path.join(kwargs['output_dir'], 'de_ssm_input.xyz')
        nreac_atoms = len(reactant.getListOfAtoms())
        nproduct_atoms = len(product.getListOfAtoms())

        with open(output, 'w') as f:
            f.write('{}\n\n{}\n{}\n\n{}\n'.format(
                nreac_atoms, reactant, nproduct_atoms, product))
        return output

    @staticmethod
    def makeCalFile(_input, filename='draw.xyz', symbol=None, geometry=None, **kwargs):
        """
        Create input file for network drawing.
        """
        output = os.path.join(kwargs['output_dir'], filename)

        if symbol is not None and geometry is not None:
            natoms = len(symbol)
            with open(output, 'w') as f:
                f.write(str(natoms))
                f.write('\n\n')
                for atom, xyz in zip(symbol, geometry):
                    f.write('{}  {}  {}  {}\n'.format(
                        atom, xyz[0], xyz[1], xyz[2]))
        else:
            natoms = len(_input.getListOfAtoms())
            with open(output, 'w') as f:
                f.write('{}\n\n{}'.format(natoms, _input))

    @staticmethod
    def makeisomerFile(add_bonds, break_bonds, **kwargs):
        """
        Create input file(add which bonds) for Single ended String Method (SSM) calculation.
        only for break 2 form 2 if more then need modify
        """
        output = os.path.join(kwargs['output_dir'], 'add_bonds.txt')

        with open(output, 'w') as f:
            if len(add_bonds) != 0:
                for i in add_bonds:
                    f.write('ADD {} {}\n'.format(i[0]+1, i[1]+1))
            if len(break_bonds) != 0:
                for i in break_bonds:
                    f.write('BREAK {} {}\n'.format(i[0]+1, i[1]+1))

    @staticmethod
    def get_mopac_opt_geometry(output, **kwargs):
        with open(output, 'r') as f:
            log = f.read().splitlines()

        iterable = reversed(range(len(log)))
        for i in iterable:
            line = log[i]
            if '                             CARTESIAN COORDINATES' in line:
                symbols, coords = [], []
                for line in log[(i+2):]:
                    if line != '':
                        data = line.split()
                        symbols.append(data[1])
                        coords.append([float(c) for c in data[2:]])
                    else:
                        return symbols, np.array(coords)

    def runXTB(self, tmpdir, config_path, target='reactant.xyz'):
        input_path = os.path.join(tmpdir, target)
        outname = '{}.xyz'.format(target.split('.')[0])
        output_path = os.path.join(tmpdir, 'xtbopt.xyz')
        new_output_path = os.path.join(tmpdir, outname)
        constraint_path = os.path.join(config_path, 'xtb_constraint.inp')

        if self.constraint == None:
            p = Popen(['xtb', input_path, '--opt', 'tight'])
            p.wait()
            os.rename(output_path, new_output_path)
        else:
            p = Popen(['xtb', '--opt', 'tight', '--input', constraint_path, input_path])
            p.wait()
            os.rename(output_path, new_output_path)

        with open(input_path, 'r') as f:
            lines = f.readlines()
        HeatofFormation = lines[1].strip().split()[1]
        return HeatofFormation