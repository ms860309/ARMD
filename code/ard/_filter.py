import os
import sys
from openbabel import openbabel as ob
from openbabel import pybel as pb
import gen3D


class FILTER(object):
    def __init__(self, reactant_file, cluster_bond_file = None, bond_file = None, fixed_atom = None):
        self.reactant_file = reactant_file
        self.cluster_bond_file = cluster_bond_file
        self.bond_file = bond_file
        self.fixed_atom = fixed_atom
        if self.fixed_atom:
            with open(self.fixed_atom, 'r') as f:
                lines = f.read()
            self.fixed_atom = eval(lines)

    def initialization(self):
        status = True
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
        status, msg = self.reactant_bonds()
        if status:
            status, msg = self.check_unreasonable_connection()
            if status:
                return 'job_success', msg
            else:
                return 'job_fail', msg
        else:
            return 'job_fail', msg

    def reactant_bonds(self):
        # extract reactant bonds
        reactant_bonds = []
        for bond in pb.ob.OBMolBondIter(self.mol.OBMol):
            sorted_bond = sorted(
                (bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1))
            bond_order = bond.GetBondOrder()
            bonds = tuple(sorted_bond) + (bond_order,)
            reactant_bonds.append(bonds)
        reactant_bonds = tuple(sorted(reactant_bonds))
        self.reactant_bonds = reactant_bonds
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
        reactant_carbon = []
        reactant_oxygen = []
        active_site_oxygen = []
        active_site_silicon = []
        active_site_metal = []  # For now only tin, tungsten and Molybdenum
        hcap = []
        for active_site_atom in self.fixed_atom:
            if self.atoms[active_site_atom] == 8:
                active_site_oxygen.append(active_site_atom)
            elif self.atoms[active_site_atom] in [42, 50, 74]:
                active_site_metal.append(active_site_atom)
            elif self.atoms[active_site_atom] == 14:
                active_site_silicon.append(active_site_atom)
            elif self.atoms[active_site_atom] == 1:
                hcap.append(active_site_atom)
        for idx, reactant_atoms in enumerate(self.atoms):
            if idx not in self.fixed_atom:
                if reactant_atoms == 8:
                    reactant_oxygen.append(idx)
                elif reactant_atoms == 6:
                    reactant_carbon.append(idx)
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


# cluster_bond = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/script/bonds.txt'
# fixed_atom = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/script/fixed_atom.txt'
# targets = """BVEMPCPWNXAEHP-UHFFFAOYSA-N_30
# BVEMPCPWNXAEHP-UHFFFAOYSA-N_72
# IJRSHDCJFOGPCF-WYGWCQTDSA-N_54
# KMHHWXAGMFWZNI-UHFFFAOYSA-M_3
# LYBKGWZNYJLYJC-XBOSCTMBSA-N_13
# NDIWZXHNLMCIMB-UHFFFAOYSA-N_12
# SXYPHYKCDMDNMD-NBINMXDDSA-N_23
# SXYPHYKCDMDNMD-WUPBIDMLSA-N_1
# """
# targets = targets.splitlines()
# for target in targets:
#     now_path = os.getcwd()
#     irc_path = os.path.join(now_path, 'irc')
#     new_path = os.path.join(irc_path, target)
#     reactant_file = os.path.join(new_path, 'irc_reactant.xyz')
#     f= FILTER(reactant_file=reactant_file, cluster_bond_file=cluster_bond, fixed_atom = fixed_atom)


# reactant_file = '/mnt/d/Lab/QMproject/AutomaticReactionDiscovery/code/ard/irc_reactant.xyz'
# f = FILTER(reactant_file=reactant_file, cluster_bond_file=cluster_bond, fixed_atom = fixed_atom)
# state, msg = f.initialization()
# print(state)
# print(msg)