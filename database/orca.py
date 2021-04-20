import os
import numpy as np


class ORCAError(Exception):
    pass


class ORCA(object):
    """
    Extract the information from ORCA output.
    The methods of this class have only been validated for ORCA 4.2.1 DFT
    calculations.
    """

    def __init__(self, outputfile=None):
        self.logfile = outputfile

        if outputfile is None:
            self.log = None
        else:
            with open(outputfile) as f:
                self.log = f.read().splitlines()
                for line in self.log:
                    if 'aborting the run' in line:
                        raise ORCAError(f'ORCA {outputfile} had an error!')
                """
                result = True
                if '                             ****ORCA TERMINATED NORMALLY****' not in self.log:
                    result = False
                if not result:
                    for line in self.log:
                        if 'aborting the run' in line:
                            raise ORCAError(f'ORCA {outputfile} had an error!')
                """

    def get_energy(self, first=False):
        # return enthalpy
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'Total Enthalpy' in line:
                return float(line.split()[-2])

    def get_entropy(self, first=False):
        # return entropy
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'Final entropy term' in line:
                return float(line.split()[-4])

    def get_free_energy(self, first=False):
        # return free_energy
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'Final Gibbs free energy' in line:
                return float(line.split()[-2])

    def get_frequencies(self):
        freqs = []
        for line in reversed(self.log):
            if 'cm**-1' in line:
                try:
                    freq = float(line.split()[1])
                except:
                    continue
                freqs.append(freq)
            elif 'VIBRATIONAL FREQUENCIES' in line:
                freqs.reverse()
                return np.array(freqs)
        else:
            raise ORCAError(f'Frequencies not found in {self.logfile}')


"""
#logfile='/mnt/d/CCEPBHYVMLWLGS-UHFFFAOYSA-N_1/TS/ts_refine.out'
logfile='/mnt/d/ts_geo.out'
try:
    q = ORCA(outputfile=logfile)
except:
    print('error')

freqs = q.get_frequencies()
nnegfreq = sum(1 for freq in freqs if freq < 0.0)
print(nnegfreq)
"""
