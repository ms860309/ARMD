import os
import sys
from openbabel import openbabel as ob
from openbabel import pybel as pb
import gen3D
import numpy as np
from statistics import mean
import props

ELEMENT_TABLE = props.ElementData()

class FILTER(object):
    def __init__(self, reactant_file, cluster_bond_file = None, fixed_atom = None):
        self.reactant_file = reactant_file
        self.cluster_bond_file = cluster_bond_file
        if fixed_atom:
            with open(self.fixed_atom, 'r') as f:
                lines = f.read()
            self.fixed_atom = eval(lines)
        
    def initialization(self):
        mol = next(pb.readfile('xyz', self.reactant_file))
        if self.cluster_bond_file:
            m = pb.ob.OBMol()
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

            with open(self.cluster_bond_file, 'r') as f:
                lines = f.read()
            cluster_bond = eval(lines)
            bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())
                    for bond in pb.ob.OBMolBondIter(mol.OBMol)]
            bonds.extend(cluster_bond)
            for bond in bonds:
                m.AddBond(bond[0], bond[1], bond[2])
            # m.ConnectTheDots()
            m.PerceiveBondOrders()
            # m.SetTotalSpinMultiplicity(1)
            m.SetTotalCharge(int(mol.charge))
            m.Center()
            m.EndModify()
            self.mol = gen3D.Molecule(m)
        else:
            self.mol = gen3D.Molecule(mol.OBMol)
        self.atoms = tuple(atom.atomicnum for atom in self.mol)

        for frag in self.mol.write('can').split()[0].split('.'):
            if '[OH]' in frag and 'Sn' not in frag:
                return 'job_fail', 'non-bonded OH group'
        return 'pass', 'pass'

    def check_feasible_rxn(self, check_mm_overlap = True, qmmm = None, qm_atoms = 23, threshold_ratio = 0.6):
        status, msg = self.initialization()
        if status == 'job_fail':
            return 'job_fail', msg

        if check_mm_overlap:
            status, msg = self.check_overlap_mm_region_v2(qmmm = qmmm, qm_atoms = qm_atoms, threshold_ratio = threshold_ratio)
        else:
            status = True
            
        if status:
            status, msg = self.check_reactant_bonds()
            if status:
                status, msg = self.check_unreasonable_connection()
                if status:
                    return 'job_success', msg
                else:
                    return 'job_fail', msg
            else:
                return 'job_fail', msg
        else:
            return 'job_fail', msg

    def check_reactant_bonds(self):
        # extract reactant bonds
        reactant_bonds = [tuple(sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)) + [bond.GetBondOrder()])
                            for bond in pb.ob.OBMolBondIter(self.mol.OBMol)]
        self.reactant_bonds = tuple(sorted(reactant_bonds))
        # check the bond order and save in a dict
        bond_type = {}
        for i in range(len(self.atoms)):
            num = 0
            for j in reactant_bonds:
                if j[0] == i or j[1] == i:
                    num += j[2]
            bond_type[i] = num
        if 0 in bond_type.values(): # The dissociated atom
            return False, 'Have dissociated atom.'
        else:
            for idx, i in enumerate(self.atoms):
                # use !=  or  >   need test
                if i == 6 and idx not in self.fixed_atom:
                    if bond_type[idx] < 3 or bond_type[idx] > 4: # remove only C=O
                        return False, 'reactant carbon bond type is invalid.({})'.format(bond_type[idx])
                # use !=  or  >   need test
                elif i == 8 and bond_type[idx] > 2 and idx not in self.fixed_atom: # Here we can't use !=2 because some time reactant O don't detect bind on Sn
                    return False, 'reactant oxygen bond type is invalid.({})'.format(bond_type[idx])
                # use !=  or  >   need test
                elif i == 14 and bond_type[idx] != 4 and idx not in self.fixed_atom:
                    return False, 'reactant silicon bond type is invalid.({})'.format(bond_type[idx])
                # While the bronsted acid already have proton on the active site, then aborted.
                elif i == 8 and bond_type[idx] > 3:
                    return False, 'oxygen have more than 4 connection.({})'.format(bond_type[idx])
            return True, 'bond type check is pass.'
    
    def check_unreasonable_connection(self):
        # Use generator is more efficient
        reactant_carbon = [idx for idx, reactant_atoms in enumerate(self.atoms) if idx not in self.fixed_atom and reactant_atoms == 6]
        reactant_oxygen = [idx for idx, reactant_atoms in enumerate(self.atoms) if idx not in self.fixed_atom and reactant_atoms == 8]
        active_site_oxygen = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] == 8]
        active_site_silicon = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] == 14]
        active_site_metal = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] in [42, 50, 74]]
        hcap = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] == 1]

        for bond in self.reactant_bonds:
            if (bond[0] in reactant_oxygen and bond[1] in active_site_silicon) or (bond[1] in reactant_oxygen and bond[0] in active_site_silicon):
                return False, 'reactant oxygen have connection with active site silicon.'
            elif (bond[0] not in self.fixed_atom and bond[1] in hcap) or (bond[1] not in self.fixed_atom and bond[0] in hcap):
                return False, 'reactant have connection with hcap.'
            elif (bond[0] in reactant_carbon and bond[1] in active_site_oxygen) or (bond[1] in reactant_carbon and bond[0] in active_site_oxygen):
                return False, 'reactant carbon have connection with active site oxygen.'
            elif (bond[0] in reactant_oxygen and bond[1] in active_site_oxygen) or (bond[1] in reactant_oxygen and bond[0] in active_site_oxygen):     
                return False, 'reactant oxygen have connection with active site oxygen.'
        return True, 'check_unreasonable_connection is pass.'

    def check_overlap_mm_region(self, qm_silicon = [], threshold = 5.4):
        # mm silicon index start from 0
        # Choose the silicon in cluster model which is in mm region
        self.mol.gen3D(self.fixed_atom, make3D=False)
        if len(self.mol.mols) == 1:
            return True, 'pass'
        else:
            nodes_1 = [mol.toNode() for mol in self.mol.mols]
            fd1 = [node.getCentroid() for node in nodes_1]
            tmps, dist2 = [], []
            for idx, i in enumerate(self.mol.mols):
                if all(idx2 not in self.fixed_atom for idx2 in i.mols_indices[idx]):
                    if any(self.atoms[idx2] == 6 or self.atoms[idx2] == 8 for idx2 in i.mols_indices[idx]):
                        tmps.append(idx)

            for tmp in tmps:
                for qm_si in qm_silicon:
                    diff2 = self.mol[qm_si].coords - fd1[tmp]
                    dist2.append(np.linalg.norm(diff2))
            # print(max(dist2))
            # print(mean(dist2))
            if max(dist2) > 6.5 and mean(dist2) > 5.4:
                # print(max(dist2))
                # print(mean(dist2))
                return False, 'Overlap with the mm region'
            else:
                return True, 'pass'

    def check_overlap_mm_region_v2(self, qmmm = None, qm_atoms = 23, threshold_ratio = 0.6):
        """
        The distance between qm atoms and mm atoms should greater than 0.6 vdw radius.
        """
        dist2 = []
        for idx1, qm_atom in enumerate(self.mol):
            if idx1 >= qm_atoms:
                continue
            for idx2, mm_atom in enumerate(qmmm):
                if idx2 < qm_atoms:
                    continue
                diff2 = np.array(qm_atom.coords) - np.array(mm_atom.coords)
                dist = np.linalg.norm(diff2)
                element = ELEMENT_TABLE.from_atomic_number(qm_atom.OBAtom.GetAtomicNum())
                vdw_rad = element.vdw_radius
                if dist < vdw_rad * threshold_ratio:
                    dist2.append(dist)
        if dist2:
            return False, 'Maybe overlap with the mm region'
        else:
            return True, 'pass'


# cluster_bond = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/script/bonds.txt'
# fixed_atom = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/script/fixed_atom.txt'
# qmmm = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/code/ard/qmmm.xyz'
# qmmm_mol = next(pb.readfile('xyz', qmmm))
# # reactant_file = os.path.join('/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/code/ard/reactions', 'UYFCJUSUJARWAH-UHFFFAOYSA-N_9/product.xyz')
# # f = FILTER(reactant_file=reactant_file, cluster_bond_file=cluster_bond, fixed_atom = fixed_atom)
# # msg = f.check_overlap_mm_region_2(qmmm = qmmm_mol, qm_atoms = 23)


# a = os.listdir('/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/code/ard/reactions')
# for i in a:
#     print('---------')
#     print(i)
#     b = os.path.join('/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/code/ard/reactions', i)
#     reactant_file = os.path.join(b, 'reactant.xyz')
#     f = FILTER(reactant_file, cluster_bond_file=cluster_bond, fixed_atom = fixed_atom)
#     status, msg = f.check_feasible_rxn(check_mm_overlap = True, qmmm = qmmm_mol, qm_atoms = 23, threshold_ratio = 0.6)
#     # state, msg = f.check_feasible_rxn(check_mm_overlap=True)
#     # print(state)
#     # print(msg)