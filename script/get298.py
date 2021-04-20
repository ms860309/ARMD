import pybel
import gen3D
from rmgpy.data.thermo import ThermoDatabase
from rmgpy import settings
import os
from os import path
import sys
sys.path.append(path.join(path.dirname(
    path.dirname(path.abspath(__file__))), 'ard'))


def get298():
    dirs = os.listdir('/mnt/d/reactions')
    for i in dirs:
        dir_path = os.path.join(os.path.join(
            '/mnt/d/reactions', i), 'product.xyz')
        mol = next(pybel.readfile('xyz', dir_path))
        mol = gen3D.Molecule(mol.OBMol)
        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'))
        thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                  'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]
        try:
            H298_prod = mol.getH298(thermo_db)
        except:
            continue

        print(i)
        print(H298_prod)


get298()
