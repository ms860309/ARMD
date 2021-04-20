class SSMError(Exception):
    pass

class SSM(object):
    def __init__(self, outputfile=None):
        self.logfile = outputfile

        if outputfile is None:
            self.log = None
        else:
            with open(outputfile) as f:
                self.log = f.read().splitlines()
    
    def check_job_status(self):
        for line in reversed(self.log):
            if 'Exiting early' in line:
                return 'Exiting early'
            elif 'Ran out of nodes' in line:
                return 'Ran out of nodes'
            elif 'termination due to dissociation' in line:
                return 'All uphill'
            elif 'total dissociation' in line:
                return 'total dissociation'
            else:
                return 'job_success'
    
    def get_ts_energy_guess(self):
        for line in reversed(self.log):
            if 'TS energy' in line:
                return float(line.split()[-1])
    
    def get_delta_e(self):
        for line in reversed(self.log):
            if 'Delta E' in line:
                return float(line.split()[-1])