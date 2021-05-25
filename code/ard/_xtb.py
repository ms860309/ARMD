# standard library imports
import os
import shutil
import time
import copy
from os import path

# third party
from subprocess import Popen, PIPE
import difflib
from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np

# local application imports
import gen3D


class XTBError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass


class XTB(object):

    def __init__(self, reactant_mol, product_mol, forcefield, constraintff_alg, form_bonds, break_bonds, logger, count, num, constraint=None, fixed_atom=None, cluster_bond = None, xtb_method = 'gfn2'):
        self.reactant_mol = reactant_mol
        self.product_mol = product_mol
        self.forcefield = forcefield
        self.constraintff_alg = constraintff_alg
        self.form_bonds = form_bonds
        self.break_bonds = break_bonds
        self.logger = logger
        self.count = count
        self.num = num
        self.constraint = constraint
        self.fixed_atom = fixed_atom
        self.cluster_bond = cluster_bond
        self.xtb_method = xtb_method

    def xtb_get_H298(self, _reactant_path, config_path):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = path.join(_reactant_path, 'tmp')
        reactant_path = path.join(tmpdir, 'reactant.xyz')
        product_path = path.join(tmpdir, 'product.xyz')

        reac_geo, prod_geo = self.genInput(self.reactant_mol, self.product_mol)

        if reac_geo == False and prod_geo == False:
            return False, False
        else:
            if path.exists(tmpdir):
                shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
            reac_geo = f"{len(str(reac_geo).splitlines())}\n\n{reac_geo}"
            prod_geo = f"{len(str(prod_geo).splitlines())}\n\n{prod_geo}"
            with open(reactant_path, 'w') as f:
                f.write(reac_geo)
            start_time = time.time()
            try:
                self.runXTB(tmpdir, config_path, constraint=self.constraint, target='reactant.xyz', method = 'gfn2')
                reactant_energy = self.getE(tmpdir, 'reactant.xyz')
            except:
                self.logger.info('xTB reactant fail')
                return False, False

            with open(product_path, 'w') as f:
                f.write(prod_geo)
            try:
                self.runXTB(tmpdir, config_path, constraint=self.constraint, target='product.xyz', method = 'gfn2')
                product_energy = self.getE(tmpdir, 'product.xyz')
            except:
                self.logger.info('xTB product fail')
                return False, False

            self.finalize(start_time, 'XTB')
            return float(reactant_energy), float(product_energy)

    def genInput(self, reactant_mol, product_mol, threshold=6.0):
        start_time = time.time()

        # Initial optimization
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField('mmff94')

        reactant_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        product_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)

        # Arrange
        try:  # Pass the more than 4 fragment situation
            arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, self.constraint, self.fixed_atom, self.cluster_bond)
            msg = arrange3D.arrangeIn3D()
            if msg != '':
                print(msg)
        except:
            self.logger.info('Here is the {} product.'.format(self.num))
            self.logger.info('Arrange fail')
            return False, False

        ff.Setup(Hatom.OBMol)
        reactant_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        ff.Setup(Hatom.OBMol)
        product_mol.gen3D(self.constraint, forcefield=self.forcefield,
                            method=self.constraintff_alg, make3D=False)
        ff.Setup(Hatom.OBMol)


        # Check reactant expected forming bond length must smaller than 4 angstrom after arrange. Default = 4
        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        # return the maximum value in array
        dist = self.check_bond_length(reactant_mol, self.form_bonds)

        if dist >= threshold:
            self.logger.info('\nHere is the {} product.'.format(self.num))
            self.logger.info('Form bonds: {}\nBreak bonds: {}\nForm bond distance: {}'.format(self.form_bonds, self.break_bonds, dist))
            self.logger.info('Form bond distance is greater than threshold.')
            self.finalize(start_time, 'arrange')
            self.logger.info('Finished {}/{}'.format(self.num, self.count))
            return False, False
        else:
            self.logger.info('\nHere is the {} product.'.format(self.num))
            self.logger.info('Structure:\n{}'.format(str(reactant_mol.toNode())))
            self.logger.info('Structure:\n{}\n'.format(str(product_mol.toNode())))
            self.logger.info('Form bonds: {}\nBreak bonds: {}\nForm bond distance: {}'.format(self.form_bonds, self.break_bonds, dist))

            product_geometry = product_mol.toNode()
            reactant_geometry = reactant_mol.toNode()
            # reactant_mol.setCoordsFromMol(reac_mol_copy)
            self.finalize(start_time, 'arrange')
            return reactant_geometry, product_geometry

    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

    @staticmethod
    def check_bond_length(product, add_bonds):
        """
        Use reactant coordinate to check if the add bonds's bond length is too long.
        Return a 'list of distance'.
        """
        atoms = tuple(atom.atomicnum for atom in product)
        coords = [atom.coords for atom in product]
        coords = [np.array(coords).reshape(len(atoms), 3)]

        dist = []
        for bond in add_bonds:
            coord_vect_1 = coords[0][bond[0]]
            coord_vect_2 = coords[0][bond[1]]
            diff = coord_vect_1 - coord_vect_2
            dist.append(np.linalg.norm(diff))

        if not dist:
            dist = [0]
        return float(max(dist))

    @staticmethod
    def getE(tmpdir, target='reactant.xyz'):
        """
        Here the energy is Eh (hartree)
        """
        input_path = path.join(tmpdir, target)
        with open(input_path, 'r') as f:
            lines = f.readlines()
        HeatofFormation = lines[1].strip().split()[1]
        return HeatofFormation

    @staticmethod
    def runXTB(tmpdir, config_path, constraint=True, target='reactant.xyz', method = 'gfn2'):
        input_path = path.join(tmpdir, target)
        outname = '{}.xyz'.format(target.split('.')[0])
        output_path = path.join(tmpdir, 'xtbopt.xyz')
        constraint_path = path.join(config_path, 'xtb_constraint.inp')

        new_output_path = path.join(tmpdir, outname)
        if constraint == None:
            if method == 'gfn2':
                p = Popen(['xtb', input_path, '--gfn', '2', '--opt', 'tight'])
            elif method == 'gfn1':
                p = Popen(['xtb', input_path, '--gfn', '1', '--opt', 'tight'])
            else:
                raise XTBError('Unsupported xtb method')
            p.wait()
            os.rename(output_path, new_output_path)
        else:
            if method == 'gfn2':
                p = Popen(['xtb', '--opt', 'tight', '--gfn', '2', '--input', constraint_path, input_path])
            elif method == 'gfn1':
                p = Popen(['xtb', '--opt', 'tight', '--gfn', '1', '--input', constraint_path, input_path])
            else:
                raise XTBError('Unsupported xtb method')
            p.wait()
            os.rename(output_path, new_output_path)
