import os
import numpy as np


class QChemError(Exception):
    pass


class QChem(object):
    """
    Extract the information from QChem output.
    The methods of this class have only been validated for Q-Chem 5.3 DFT
    calculations.
    """

    def __init__(self, outputfile=None):
        self.logfile = outputfile

        if outputfile is None:
            self.log = None
        else:
            with open(outputfile) as f:
                self.log = f.read().splitlines()
                result = True
                restart = True
                special_case = True
                if ' **  OPTIMIZATION CONVERGED  **' not in self.log:
                    result = False
                if ' **  MAXIMUM OPTIMIZATION CYCLES REACHED  **' not in self.log:
                    restart = False
                if ' ** UNABLE TO DETERMINE Lamda IN FormD **' not in self.log:
                    special_case = True

                if not result and not restart and not special_case:
                    for line in self.log:
                        if 'fatal error' in line:
                            raise QChemError(f'Q-Chem job {outputfile} had an error!')

    def get_energy(self, first=False):
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'total energy' in line:  # Double hybrid methods
                return float(line.split()[-2])
            elif 'energy in the final basis set' in line:  # Other DFT methods
                return float(line.split()[-1])
            elif 'The QM part of the Energy is' in line: # QMMM
                return float(line.split()[-1])
        else:
            raise QChemError(f'Energy not found in {self.logfile}')

    def get_geometry(self, first=False):
        if first:
            iterable = range(len(self.log))
        else:
            iterable = reversed(range(len(self.log)))
        for i in iterable:
            line = self.log[i]
            if 'Standard Nuclear Orientation' in line:
                symbols, coords = [], []
                for line in self.log[(i+3):]:
                    if '----------' not in line:
                        data = line.split()
                        symbols.append(data[1])
                        coords.append([float(c) for c in data[2:]])
                    else:
                        return symbols, np.array(coords)
        else:
            raise QChemError(f'Geometry not found in {self.logfile}')

    def get_frequencies(self):
        freqs = []
        for line in reversed(self.log):
            if 'Frequency' in line:
                freqs.extend([float(f) for f in reversed(line.split()[1:])])
            elif 'VIBRATIONAL ANALYSIS' in line:
                freqs.reverse()
                return np.array(freqs)
        else:
            raise QChemError(f'Frequencies not found in {self.logfile}')

    def get_normal_modes(self):
        modes = []
        for i in reversed(range(len(self.log))):
            line = self.log[i]
            if 'Raman Active' in line:
                mode1, mode2, mode3 = [], [], []
                for line in self.log[(i+2):]:
                    if 'TransDip' not in line:
                        vals = line.split()[1:]
                        mode1.append([float(v) for v in vals[:3]])
                        mode2.append([float(v) for v in vals[3:6]])
                        mode3.append([float(v) for v in vals[6:]])
                    else:
                        modes.extend(
                            [np.array(mode3), np.array(mode2), np.array(mode1)])
                        break
            elif 'VIBRATIONAL ANALYSIS' in line:
                modes.reverse()
                return modes
        else:
            raise QChemError(f'Normal modes not found in {self.logfile}')

    def get_zpe(self):
        for line in reversed(self.log):
            if 'Zero point vibrational energy' in line:
                return float(line.split()[-2]) / 627.5095  # Convert to Hartree
        else:
            raise QChemError(f'ZPE not found in {self.logfile}')

    def create_geo_file(self, file_path):
        symobol, geometry = self.get_geometry()
        natoms = len(symobol)
        with open(file_path, 'w') as f:
            f.write(str(natoms))
            f.write('\n\n')
            for atom, xyz in zip(symobol, geometry):
                f.write('{}  {}  {}  {}\n'.format(atom, xyz[0], xyz[1], xyz[2]))

    def get_opt_cycle(self):
        count = 0
        for i in self.log:
            if i == '   Searching for a Minimum':
                count += 1
        return count
    
    def create_qmmm_geomtry(self, file_path):
        qm_atoms = []
        for i, text in enumerate(self.log):
            if text.upper().startswith('$QM_ATOMS'):
                break
        for qm_atom in self.log[i+1:]:
            if '$end' in qm_atom:
                break
            else:
                try:
                    qm_atoms.append(int(qm_atom))
                except:
                    continue
        symobol, geometry = self.get_geometry()
        natoms = len(qm_atoms)
        with open(file_path, 'w') as f:
            f.write(str(natoms))
            f.write('\n\n')
            for atom, xyz in zip(symobol[:natoms], geometry[:natoms]):
                f.write('{}  {}  {}  {}\n'.format(atom, xyz[0], xyz[1], xyz[2]))


# logfile='/mnt/d/Lab/QMproject/AutomatedReactionMechanismDiscovery/database/qmmm_freq_ts.out'
# try:
#     q = QChem(outputfile=logfile)
# except:
#     print('error')
