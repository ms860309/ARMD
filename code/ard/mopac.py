# standard library imports
import os
from os import path
import shutil
import time

# third party
from subprocess import Popen, PIPE
import difflib
from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np

# local application imports
from node import Node
import gen3D
import util


class MopacError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass


class Mopac(object):

    def __init__(self, reactant_mol, product_mol, mopac_method, forcefield, constraintff_alg, form_bonds, break_bonds, logger, count, num, constraint=None, fixed_atom=None, cluster_bond = None):
        self.reactant_mol = reactant_mol
        self.product_mol = product_mol
        self.mopac_method = mopac_method
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

    def mopac_get_H298(self, tmp_path, charge=0, multiplicity='SINGLET'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = path.join(tmp_path, 'tmp')
        reactant_path = path.join(tmpdir, 'reactant.mop')
        product_path = path.join(tmpdir, 'product.mop')

        reac_geo, prod_geo = self.genInput(self.reactant_mol, self.product_mol)

        if reac_geo == False and prod_geo == False:
            self.logger.info('Mopac fail')
            return False, False
        else:
            if path.exists(tmpdir):
                shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)
            with open(reactant_path, 'w') as f:
                f.write("NOSYM CHARGE={} {} {}\n\n".format(charge, multiplicity, self.mopac_method))
                f.write("\n{}".format(reac_geo))
            start_time = time.time()
            try:
                self.runMopac(tmpdir, 'reactant.mop')
                reactant = self.getHeatofFormation(tmpdir, 'reactant.out')
            except:
                self.logger.info('Mopac reactant fail')
                return False, False

            with open(product_path, 'w') as f:
                f.write("NOSYM CHARGE={} {} {}\n\n".format(charge, multiplicity, self.mopac_method))
                f.write("\n{}".format(prod_geo))
            try:
                self.runMopac(tmpdir, 'product.mop')
                product = self.getHeatofFormation(tmpdir, 'product.out')
            except:
                self.logger.info('Mopac product fail')
                return False, False
            self.finalize(start_time, 'mopac')
            return float(reactant), float(product)

    def genInput(self, reactant_mol, product_mol, threshold=5.0):
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
            arrange3D = gen3D.Arrange3D(
                reactant_mol, product_mol, self.constraint, self.fixed_atom)
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
            
            reactant_geometry = self.gen_geo_inp(reactant_mol, constraint=self.constraint)
            product_geometry = self.gen_geo_inp(product_mol, constraint=self.constraint)

            self.finalize(start_time, 'arrange')
            return reactant_geometry, product_geometry

    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

    @staticmethod
    def gen_geo_inp(mol_object, constraint=None):
        geometry = str(mol_object.toNode()).splitlines()
        product_geometry = []
        for idx, i in enumerate(geometry):
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            if not constraint:
                l = " 1 ".join(k)
            else:
                if idx in constraint:
                    l = " 0 ".join(k)
                else:
                    l = " 1 ".join(k)
            out = atom + l
            product_geometry.append(out)
        product_geometry = "\n".join(product_geometry)
        return product_geometry

    @staticmethod
    def check_bond_length(product, add_bonds):
        """
        Use reactant coordinate to check if the add bonds's bond length is too long.
        Return a 'list of distance'.
        """
        coords = [atom.coords for atom in product]
        atoms = tuple(atom.atomicnum for atom in product)
        coords = [np.array(coords).reshape(len(atoms), 3)]

        dist = []
        for bond in add_bonds:
            coord_vect_1 = coords[0][bond[0]]
            coord_vect_2 = coords[0][bond[1]]
            diff = coord_vect_1 - coord_vect_2
            dist.append(np.linalg.norm(diff))

        if dist == []:
            dist = [0]
        return float(max(dist))

    @staticmethod
    def getHeatofFormation(tmpdir, target='reactant.out'):
        """
        if Error return False, which HF may be 0.0
        """
        input_path = path.join(tmpdir, target)
        with open(input_path, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.strip().startswith('FINAL HEAT OF FORMATION'):
                break
        string = lines[idx].split()
        if string[0] == 'FINAL':
            HeatofFormation = string[5]
        else:
            HeatofFormation = False
        return HeatofFormation

    @staticmethod
    def runMopac(tmpdir, target='reactant.mop'):
        input_path = path.join(tmpdir, target)
        p = Popen(['mopac', input_path])
        p.wait()
