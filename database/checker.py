# connect to database
# search entries with status "job_launched" or "job_running"
# check the job pbs-status, e.g., qstat -f "job_id"
# 1. Invalid job id specified: means the job is off queue
#    1.1 check job folder if input.log file is there
#			Yes: job already finished
#			No:  job is aborted before even calculation start: job_aborted
#    1.2 if finished, check "Normal termination" key words should appear twice
#			Yes: job is converged
#			No:  job fails convergence (possibly need more wall-time): job_failed_convergence
#    1.3 if converged, check isomorphism between input molecule and output molecule
#			Yes: job is success: job_success
#			No:  job fails isomorphism (possibly need tweak on initial structure): job_failed_isomorphism
# 2. output string starts with e.g, "JobId=5037088"
#    2.1 check JobState
#		RUNNING: running: job_running
#		else:   still job_launched

from connect import db
import subprocess
import os
from os import path
import sys
import shutil
from openbabel import pybel
from openbabel import openbabel as ob
from qchem import QChem
from orca import ORCA
from ssm import SSM
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'code/ard'))
from _filter import FILTER
from typing import Union

class CheckError(Exception):
    pass

"""
Global utils
"""

def select_targets(qm_collection:object, job_name:str) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched, job_running or job_queueing
    Returns a list of target
    """
    keyname = '{}_status'.format(job_name)
    query = {keyname:
             {"$in":
              ["job_launched", "job_running", "job_queueing"]
              }
             }
    targets = list(qm_collection.find(query))
    return targets

def check_job_status(job_id:str) -> str:
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """

    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if "Unknown Job Id" in stderr.decode():
        return "off_queue"

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return "job_running"
    elif stdout[idx+2] == 'Q':
        return 'job_queueing'
    else:
        return "job_launched"

def xyz_to_pyMol(xyz:str, cluster_bond_path:str=None) -> object:
    mol = next(pybel.readfile('xyz', xyz))
    if cluster_bond_path:
        m = pybel.ob.OBMol()
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

        with open(cluster_bond_path, 'r') as f:
            lines = f.read()
        cluster_bond = eval(lines)
        bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())
                 for bond in pybel.ob.OBMolBondIter(mol.OBMol)]
        bonds.extend(cluster_bond)
        for bond in bonds:
            m.AddBond(bond[0], bond[1], bond[2])
        # m.ConnectTheDots()
        m.PerceiveBondOrders()
        # m.SetTotalSpinMultiplicity(1)
        m.SetTotalCharge(int(mol.charge))
        m.Center()
        m.EndModify()

        pybelmol = pybel.Molecule(m)
        return pybelmol
    else:
        return mol

"""
SSM check.
"""

def check_ssm_content(target_path:str, thershold:float = 200.0) -> str:
    """
    Check ssm log file. Will have three status
    1. job success
    2. job fail (All uphill, totol dissociation, etc.  see ssm.py)
    3. ts node high energy
    """
    status_path = path.join(target_path, 'SSM/status.log')
    ts_node_path = path.join(target_path, 'SSM/TSnode.xyz')
    try:
        q = SSM(outputfile=status_path)
        job_status = q.check_job_status()
        ts_energy_guess = q.get_ts_energy_guess()
        # delta_e_guess = q.get_delta_e()
        if ts_energy_guess > thershold:
            return "ts_guess high energy"
        elif os.path.exists(ts_node_path):
            return job_status
        else:
            return "job_fail"
    except:
        return "job_fail"

def check_ssm_jobs(qm_collection:object, refine:bool=False, thershold:float = 200.0):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='ssm')
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ssm_jobid']
        # 2. check the job pbs status
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_ssm_content(target['path'], thershold = thershold)
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ssm_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                if refine:
                    update_field = {
                        'ssm_status': new_status, "ts_refine_status": "job_unrun"
                        }
                else:
                    update_field = {
                        'ssm_status': new_status, "ts_status": "job_unrun"
                        }
            else:
                update_field = {
                    'ssm_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
TS refine check
"""

def check_ts_refine_content(target_path:str, threshold:float=-50.0) -> str:
    """
    Check ts refine log file. Will have three status
    1. job success
    2. job crash
    3. imaginary frequency do not equal to 1 or greater than -50 cm^-1.
    """
    ts_dir_path = path.join(target_path, 'TS')
    ts_refine_out_path = path.join(ts_dir_path, 'ts_refine.out')

    if os.path.exists(ts_refine_out_path):
        try:
            q = ORCA(outputfile=ts_refine_out_path)
            freqs = q.get_frequencies()
            nnegfreq = sum(1 for freq in freqs if freq < 0.0)
            if nnegfreq > 1:
                return "Have more than one imaginary frequency"
            elif nnegfreq == 0:
                return "All positive frequency"
            else:
                if min(freqs) < threshold:
                    #ts_energy = q.get_free_energy()
                    #ts_energy = q.get_energy()
                    return "job_success"
                else:
                    return f"Imaginary frequency greater than {threshold} cm-1"
        except:
            return "job_fail"
    else:
        return "job_fail"

def check_ts_refine_jobs(qm_collection:object, threshold:float=-50.0):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='ts_refine')
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ts_refine_jobid']
        # 2. check the job pbs status
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            new_status = check_ts_refine_content(target['path'], threshold = threshold)
        orig_status = target['ts_refine_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                   'ts_status':'job_unrun','ts_refine_status': new_status
                    }
            else:
                update_field = {
                    'ts_refine_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
TS check.
"""

def check_ts_content(target_path:str, threshold:float = -50.0) -> Union[str, float]:
    level_of_theory = 'QCHEM'
    ts_dir_path = path.join(target_path, 'TS')
    qchem_ts_out_path = path.join(ts_dir_path, 'ts.out')
    orca_ts_out_path = path.join(ts_dir_path, 'ts_geo.out')
    # Which means the ts is running with orca
    if not os.path.exists(qchem_ts_out_path):
        level_of_theory = 'ORCA'
    ts_geo_path = path.join(ts_dir_path, 'ts_geo.xyz')

    if level_of_theory == 'QCHEM':
        try:
            q = QChem(outputfile=qchem_ts_out_path)
            freqs = q.get_frequencies()
            nnegfreq = sum(1 for freq in freqs if freq < 0.0)
            if nnegfreq > 1:
                return "Have more than one imaginary frequency", 0.0
            elif nnegfreq == 0:
                return "All positive frequency", 0.0
            else:
                if min(freqs) < threshold:
                    energy = q.get_energy()
                    zpe = q.get_zpe()
                    ts_energy = energy + zpe
                    q.create_geo_file(ts_geo_path)
                    return "job_success", float(ts_energy)
                else:
                    return f"Imaginary frequency greater than {threshold} cm-1", 0.0
        except:
            return "job_fail", 0.0
    elif level_of_theory == 'ORCA':
        if os.path.exists(ts_geo_path):
            try:
                q = ORCA(outputfile=orca_ts_out_path)
                freqs = q.get_frequencies()
                nnegfreq = sum(1 for freq in freqs if freq < 0.0)
                if nnegfreq > 1:
                    return "Have more than one imaginary frequency", 0.0
                elif nnegfreq == 0:
                    return "All positive frequency", 0.0
                else:
                    if min(freqs) < threshold:
                        #ts_energy = q.get_free_energy()
                        ts_energy = q.get_energy()
                        return "job_success", float(ts_energy)
                    else:
                        return f"Imaginary frequency greater than {threshold} cm-1", 0.0
            except:
                return "job_fail", 0.0
        else:
            return "job_fail", 0.0

def check_ts_jobs(qm_collection:object, threshold:float = -50.0):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='ts')
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ts_jobid']
        # 2. check the job pbs status
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, ts_energy = check_ts_content(target['path'], threshold = threshold)
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ts_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                if target['use_irc'] == '0':
                    update_field = {
                        'ts_status': new_status, 'ts_energy': ts_energy, 'energy_status': 'job_unrun'
                        }
                else:
                    update_field = {
                        'ts_status': new_status, 'ts_energy': ts_energy, 'irc_status': 'job_unrun'
                        }
            else:
                update_field = {
                    'ts_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
IRC check.
"""

def check_irc_content(target_path:str) -> str:
    irc_dir_path = os.path.join(target_path, 'IRC')
    first_output = os.path.join(irc_dir_path, 'finished_first.xyz')
    last_output = os.path.join(irc_dir_path, 'finished_last.xyz')
    forward_end_output = os.path.join(irc_dir_path, 'forward_end_opt.xyz')
    backward_end_output = os.path.join(irc_dir_path, 'backward_end_opt.xyz')
    if os.path.exists(forward_end_output) and os.path.exists(backward_end_output):
        return 'job_success'
    # The last two is to make sure not in the crash
    elif os.path.exists(first_output) and os.path.exists(last_output) and not os.path.exists(forward_end_output) and not os.path.exists(backward_end_output):
        return 'job_success'
    else:
        return 'job_fail'

def check_irc_jobs(qm_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='irc')
    # 2. check the job pbs status
    for target in targets:
        job_id = target['irc_jobid']
        # 2. check the job pbs status
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_content(target['path'])
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['irc_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'irc_status': new_status, 'irc_opt_status': 'job_unrun'
                    }
            else:
                update_field = {
                    'irc_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)


"""
IRC  success check.
This is to check whether irc forward and backward direction equal to expectation.
If one direction equal to reactant inchi key then it's an intended reaction.
"""

def select_irc_equal_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {"irc_equal":
             {"$in":
              ["waiting for checking"]
              }
             }
    targets = list(qm_collection.find(query))
    return targets

def check_irc_equal(qm_collection:object, cluster_bond_path:str=None, fixed_atom_path:str=None):
    targets = select_irc_equal_target(qm_collection)
    acceptable_condition = ['forward equal to reactant',
                            'backward equal to reactant']
    # special_condition = ['forward equal to reverse', 'unintended']

    for target in targets:
        new_status, forward, backward = check_irc_equal_status(target, cluster_bond_path=cluster_bond_path, fixed_atom_path = fixed_atom_path)
        orig_status = target['irc_equal']
        if orig_status != new_status:
            if new_status in acceptable_condition:
                update_field = {
                    'irc_equal': new_status, 'barrier': 'need check', 'insert_reaction': 'need insert',
                    'reactant_inchi_key': forward.write('inchiKey').strip(), 'product_inchi_key': backward.write('inchiKey').strip(),
                    'reactant_smiles': forward.write('can').split()[0], 'product_smiles': backward.write('can').split()[0]
                }
            else:
                update_field = {
                    'irc_equal': new_status,
                    'reactant_inchi_key': forward.write('inchiKey').strip(), 'product_inchi_key': backward.write('inchiKey').strip(),
                    'reactant_smiles': forward.write('can').split()[0], 'product_smiles': backward.write('can').split()[0]
                }
            qm_collection.update_one(target, {"$set": update_field}, True)

def check_irc_equal_status(target:object, cluster_bond_path:str=None, fixed_atom_path:str=None) -> Union[str, object, object]:
    irc_path = path.join(target['path'], 'IRC/')
    forward_end_output = os.path.join(irc_path, 'irc_forward.xyz')
    backward_end_output = os.path.join(irc_path, 'irc_backward.xyz')
    irc_reactant_path = path.join(target['path'], 'irc_reactant.xyz')  # for next generation
    pyMol_3 = xyz_to_pyMol(forward_end_output, cluster_bond_path=cluster_bond_path)
    pyMol_4 = xyz_to_pyMol(backward_end_output, cluster_bond_path=cluster_bond_path)

    reactant_inchi_key = target['reactant_inchi_key']
    if pyMol_3.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward equal to reverse', pyMol_3, pyMol_4
    elif pyMol_3.write('inchiKey').strip() == reactant_inchi_key:
        f = FILTER(reactant_file=backward_end_output, cluster_bond_file=cluster_bond_path, fixed_atom = fixed_atom_path)
        status, msg = f.initialization()
        f2 = FILTER(reactant_file=forward_end_output, cluster_bond_file=cluster_bond_path, fixed_atom = fixed_atom_path)
        status2, msg2 = f2.initialization()
        if status == 'job_success' and status2 == 'job_success':
            shutil.copyfile(backward_end_output, irc_reactant_path)
            return 'forward equal to reactant', pyMol_3, pyMol_4
        else:
            return msg, pyMol_3, pyMol_4
    elif pyMol_4.write('inchiKey').strip() == reactant_inchi_key:
        f = FILTER(reactant_file=forward_end_output, cluster_bond_file=cluster_bond_path, fixed_atom = fixed_atom_path)
        status, msg = f.initialization()
        f2 = FILTER(reactant_file=backward_end_output, cluster_bond_file=cluster_bond_path, fixed_atom = fixed_atom_path)
        status2, msg2 = f2.initialization()
        if status == 'job_success' and status2 == 'job_success':
            shutil.copyfile(forward_end_output, irc_reactant_path)
            return 'backward equal to reactant', pyMol_4, pyMol_3
        else:
            return msg, pyMol_4, pyMol_3
    else:
        return 'unintended', pyMol_3, pyMol_4

"""
IRC opt check.
This check is to make irc reactant as a new reactant to run ard next generation
"""

def select_irc_opt_finished_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$and':
             [
                 {"irc_forward_opt_status":
                  {"$in":
                   ['job_success']}},
                 {'irc_backward_opt_status':
                  {'$in':
                   ['job_success']}}
             ]
             }
    targets = list(qm_collection.find(query))
    return targets


def check_irc_opt_content(dir_path:str, level_of_theory:str='ORCA', direction:str='forward') -> Union[str, float]:
    irc_path = path.join(dir_path, "IRC")
    if level_of_theory == 'QCHEM':
        if direction == 'forward':
            output_path = path.join(irc_path, 'irc_forward.out')
        else:
            output_path = path.join(irc_path, 'irc_backward.out')
        try:
            q = QChem(outputfile=output_path)
            freqs = q.get_frequencies()
            nnegfreq = sum(1 for freq in freqs if freq < 0.0)
            if nnegfreq > 0:
                return 'Have negative frequency', 0.0
            else:
                energy = q.get_energy()
                zpe = q.get_zpe()
                energy += zpe
                return 'job_success', energy
        except:
            return 'job_fail', 0.0
    elif level_of_theory == 'ORCA':
        if direction == 'forward':
            output_path = path.join(irc_path, 'irc_forward.out')
        else:
            output_path = path.join(irc_path, 'irc_backward.out')
        try:
            q = ORCA(outputfile=output_path)
            freqs = q.get_frequencies()
            nnegfreq = sum(1 for freq in freqs if freq < 0.0)
            if nnegfreq > 0:
                return 'Have negative frequency', 0.0
            else:
                #energy = q.get_free_energy()
                energy = q.get_energy()
                return 'job_success', float(energy)
        except:
            return 'job_fail', 0.0

def check_irc_opt_jobs(qm_collection:object, level_of_theory:str='ORCA'):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='irc_forward_opt')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['irc_forward_opt_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_irc_opt_content(target['path'], level_of_theory=level_of_theory.upper(), direction='forward')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['irc_forward_opt_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'irc_forward_opt_status': new_status, 'irc_forward_opt_energy': energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'irc_forward_opt_status': new_status
                    }
            else:
                update_field = {
                    'irc_forward_opt_status': new_status, 'irc_forward_opt_energy': energy
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='irc_backward_opt')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['irc_backward_opt_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_irc_opt_content(target['path'], level_of_theory=level_of_theory.upper(), direction='backward')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['irc_backward_opt_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'irc_backward_opt_status': new_status, 'irc_backward_opt_energy': energy,
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'irc_backward_opt_status': new_status
                    }
            else:
                update_field = {
                    'irc_backward_opt_status': new_status, 'irc_backrward_opt_energy': energy
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    finished_targets = select_irc_opt_finished_target(qm_collection)
    for target in finished_targets:
        update_field = {
            'irc_opt_status': 'job_success', 'irc_equal': 'waiting for checking'
            }
        qm_collection.update_one(target, {"$unset": {'irc_backward_opt_status': '', 'irc_forward_opt_status': '',
                                                     'irc_backward_opt_jobid': '', 'irc_forward_opt_jobid': ''}, "$set": update_field}, True)

"""
IRC opt side fail check.
If one side irc opt fail then delete the other side (while the other side still running).
"""

def select_irc_opt_side_fail_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$or':[{'$and':
            [
                {"irc_forward_opt_status":
                {"$in":
                    ['job_fail', 'Have negative frequency']}},
                {'irc_backward_opt_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            },
            {'$and':
            [
                {"irc_backward_opt_status":
                {"$in":
                    ['job_fail', 'Have negative frequency']}},
                {'irc_forward_opt_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            }]}
    targets = list(qm_collection.find(query))
    return targets

def check_irc_opt_side_fail_jobs(qm_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_irc_opt_side_fail_target(qm_collection)
    # 2. check the job pbs_status
    for target in targets:
        if target['irc_forward_opt_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['irc_backward_opt_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'irc_forward_opt_status': 'the other side fail, so this side be deleted'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['irc_backward_opt_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['irc_forward_opt_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'irc_backward_opt_status': 'the other side fail, so this side be deleted'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
After irc check, insert the reaction into reaction collection.
Here we select the lowest activation energy one.  #TO DO
If we insert to reaction collection, next generation is ready to run. (In the other words, create ard_status: job_unrun)  <-- the lowest activation energy one
BTW, when insert reaction, we may consider the same reaction had been generated by early generation.
# To make sure the lowest barrier add to reaction collection.
# So we need to know all of the job are finished.
# The basic thought is that there is not any unrun job except ard job.
# def create ard job  #TO DO
"""

def select_insert_reaction_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$and':
             [
                 {"insert_reaction":
                  {"$in":
                   ['need insert']}},
                 {'barrier':
                  {'$nin':
                   ['need check']}}]
             }
    targets = list(qm_collection.find(query))
    return targets


def insert_reaction(qm_collection:object, reactions_collection:object):
    targets = select_insert_reaction_target(qm_collection)

    acceptable_condition = ['forward equal to reactant',
                            'backward equal to reactant']

    # new one not mean the lowest barrier (so the lowest may in duplicate)
    for target in targets:
        reactant_inchi_key = target['reactant_inchi_key']
        product_inchi_key = target['product_inchi_key']
        reactant_smiles = target['reactant_smiles']
        product_smiles = target['product_smiles']
        path = target['path']
        generations = target['generations']
        barrier = target['barrier']
        check1 = {'reaction': [reactant_inchi_key, product_inchi_key]}
        checker1 = list(reactions_collection.find(check1))
        irc_equal = target['irc_equal']

        if irc_equal in acceptable_condition:
            if target['irc_equal'] == 'forward equal to reactant':
                delta_H = (target['irc_backward_opt_energy'] - target['irc_forward_opt_energy']) * 627.5095
            elif target['irc_equal'] == 'backward equal to reactant':
                delta_H = (target['irc_forward_opt_energy'] - target['irc_backward_opt_energy']) * 627.5095
            if len(checker1) == 0:
                reactions_collection.insert_one({
                    'reaction': [reactant_inchi_key, product_inchi_key],
                    'reactant_inchi_key': reactant_inchi_key,
                    'product_inchi_key': product_inchi_key,
                    'reactant_smiles': reactant_smiles,
                    'product_smiles': product_smiles,
                    'path': path,
                    'generations': generations,
                    'unique': 'new one',
                    'barrier': barrier,
                    'irc_equal': irc_equal,
                    'irc_opt_status': target['irc_opt_status'],
                    'delta_H': delta_H})
            else:
                reactions_collection.insert_one({
                    'reaction': [reactant_inchi_key, product_inchi_key],
                    'reactant_inchi_key': reactant_inchi_key,
                    'product_inchi_key': product_inchi_key,
                    'reactant_smiles': reactant_smiles,
                    'product_smiles': product_smiles,
                    'path': path,
                    'generations': generations,
                    'unique': 'duplicate',
                    'barrier': barrier,
                    'irc_equal': irc_equal,
                    'irc_opt_status': target['irc_opt_status'],
                    'delta_H': delta_H})
        else:
            if len(checker1) == 0:
                reactions_collection.insert_one({
                    'reaction': [reactant_inchi_key, product_inchi_key],
                    'reactant_inchi_key': reactant_inchi_key,
                    'product_inchi_key': product_inchi_key,
                    'reactant_smiles': reactant_smiles,
                    'product_smiles': product_smiles,
                    'path': path,
                    'generations': generations,
                    'barrier': barrier,
                    'ard_status': 'already insert to qm',
                    'irc_equal': irc_equal})
            else:
                reactions_collection.insert_one({
                    'reaction': [reactant_inchi_key, product_inchi_key],
                    'reactant_inchi_key': reactant_inchi_key,
                    'product_inchi_key': product_inchi_key,
                    'reactant_smiles': reactant_smiles,
                    'product_smiles': product_smiles,
                    'path': path,
                    'generations': generations,
                    'unique': 'duplicate',
                    'barrier': barrier,
                    'ard_status': 'already insert to qm',
                    'irc_equal': irc_equal})
        qm_collection.update_one(target, {"$set": {'insert_reaction': "already insert to qm"}}, True)

"""
ARD check unrun
"""

def insert_ard(qm_collection:object, reactions_collection:object, statistics_collection:object, config_collection:object, barrier_threshold:float=60.0):
    use_irc = list(config_collection.find({'generations':1}))[0]['use_irc']
    ard_query = {"ard_status":
                 {"$in":
                  ["job_unrun", "job_launched", "job_running", "job_queueing"]
                  }
                 }
    energy_query = {"energy_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                     }
                    }
    ssm_query = {"ssm_status":
                 {"$in":
                  ["job_unrun", "job_launched", "job_running", "job_queueing"]
                  }
                 }
    ts_refine_query = {"ts_refine_status":
                       {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                        }
                       }
    ts_query = {"ts_status":
                {"$in":
                 ["job_unrun", "job_launched", "job_running", "job_queueing"]
                 }
                }
    irc_query = {"irc_status":
                 {"$in":
                  ["job_unrun", "job_launched", "job_running", "job_queueing"]
                  }
                 }
    irc_opt_query = {'$or':
                     [
                         {"irc_forward_opt_status":
                          {"$in":
                           ["job_unrun", "job_launched", "job_running", "job_queueing"]}},
                         {'irc_forward_opt_status':
                             {'$in':
                              ["job_unrun", "job_launched", "job_running", "job_queueing"]}},
                         {'irc_opt_status':
                             {'$in':
                              ["job_unrun"]}}
                     ]
                     }
    irc_equal_query = {"irc_equal":
                       {"$in":
                        ["waiting for checking"]
                        }
                       }
    insert_reaction_query = {"insert_reaction":
                             {"$in":
                              ['need insert']}}

    if use_irc == '0':
        not_finished_number = len(list(qm_collection.find({'$or':
                                                [energy_query, ssm_query, ts_query, insert_reaction_query, ts_refine_query, ard_query]
                                                })))
    else:
        not_finished_number = len(list(qm_collection.find({'$or':
                                                [energy_query, ssm_query, ts_query, irc_query, irc_opt_query, irc_equal_query, insert_reaction_query, ts_refine_query, ard_query]
                                                })))

    ard_had_add_number = qm_collection.count_documents({})
    ard_should_add_number = 0
    should_adds = list(statistics_collection.find({}))
    for i in should_adds:
        ard_should_add_number += i['add how many products']

    if int(not_finished_number) == 0 and int(ard_had_add_number) == int(ard_should_add_number):
        insert_qmmm(qm_collection, reactions_collection)

        acceptable_condition = ['forward equal to reactant',
                                'backward equal to reactant']

        finished_reactant_list = []
        finished_reactant_smiles_part_list = []
        for i in should_adds:
            finished_reactant_list.append(i['reactant_inchi_key'])
            reactant_smiles = i['reactant_smiles'].split('.')
            if len(reactant_smiles) > 1:
                reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs and 'C' in rs])
            finished_reactant_smiles_part_list.append(reactant_part_smiles)

        reactions = list(reactions_collection.aggregate([{
            '$match': {
                'irc_opt_status': 'job_success',
                'irc_equal': {'$in': acceptable_condition},
                'ard_status': {'$nin': ['already insert to qm']},
                'product_inchi_key':{'$nin': finished_reactant_list},
                'barrier': {'$lte': barrier_threshold, '$gte': -100.0},
                # Filter the energy fail
                'delta_H': {'$lt': 200, '$gte': -10000.0} #-10000 is to filter the failed jobs
            }}, {
            '$group': {
                '_id': "$reaction",
                'barrier': {'$min': "$barrier"}
            }}
        ]))
        # If the same product but different active site, then choose the lowest barrier one.
        tmp = {}
        for i in reactions:
            dir_path = list(reactions_collection.find({'reaction': i['_id'], 'barrier': i['barrier']}))[0]
            ard_qm_target = list(qm_collection.find({'path': dir_path['path']}))[0]
            # Prevent reactant equal to product but with different active site (maybe proton at the different oxygen)
            same = False
            reactant_smiles = ard_qm_target['reactant_smiles'].split('.')
            product_smiles = ard_qm_target['product_smiles'].split('.')
            if len(reactant_smiles) > 1:
                reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs and 'C' in rs])
            else:
                reactant_part_smiles = set(reactant_smiles)

            if len(product_smiles) > 1:
                product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps and 'C' in ps])
            else:
                product_part_smiles = set(product_smiles)

            if reactant_part_smiles == product_part_smiles:
                same = True
            # Prevent the product already be a reactant before
            if not same:
                for finished_reactant_part in finished_reactant_smiles_part_list:
                    if product_part_smiles == finished_reactant_part:
                        same = True
                        break
            if not same and ard_qm_target['product_smiles'] not in tmp:
                tmp[ard_qm_target['product_smiles']] = [i['barrier'], ard_qm_target, dir_path]
            elif not same and tmp[ard_qm_target['product_smiles']][0] > i['barrier']:
                tmp[ard_qm_target['product_smiles']] = [i['barrier'], ard_qm_target, dir_path]

        for i in tmp:
            update_field_for_qm_target = {'ard_status': 'job_unrun'}
            update_field_for_reaction_target = {'ard_status': 'already insert to qm'}
            qm_collection.update_one(tmp[i][1], {"$set": update_field_for_qm_target}, True)
            reactions_collection.update_one(tmp[i][2], {"$set": update_field_for_reaction_target}, True)

"""
QMMM
"""

def insert_qmmm(qm_collection:object, reactions_collection:object):

    reactions = list(reactions_collection.aggregate([{
                '$group': {
                    '_id': "$reaction",
                    'barrier': {'$min': "$barrier"}
                }}
            ]))

    for i in reactions:
        dir_path = list(reactions_collection.find({'reaction': i['_id'], 'barrier': i['barrier']}))[0]
        try:
            qmmm = dir_path['qmmm']
            if qmmm == 'Already insert':
                continue
        except:
            pass
        ard_qm_target = list(qm_collection.find({'path': dir_path['path']}))[0]
        # Prevent reactant equal to product but with different active site (maybe proton at the different oxygen)
        reactant_smiles = ard_qm_target['reactant_smiles'].split('.')
        product_smiles = ard_qm_target['product_smiles'].split('.')
        if len(reactant_smiles) > 1:
            reactant_part_smiles = set([rs for rs in reactant_smiles if 'Sn' not in rs and 'C' in rs])
        else:
            reactant_part_smiles = set(reactant_smiles)

        if len(product_smiles) > 1:
            product_part_smiles = set([ps for ps in product_smiles if 'Sn' not in ps and 'C' in ps])
        else:
            product_part_smiles = set(product_smiles)

        if reactant_part_smiles == product_part_smiles:
            continue

        qm_collection.update_one(ard_qm_target, {"$set": {"qmmm_opt_status": "job_unrun", "qmmm_freq_ts_status": "job_unrun", 'qmmm_freq_opt_reactant_restart_times':0, 'qmmm_freq_opt_product_restart_times':0, 'qmmm_freq_ts_restart_times':0}}, True)
        reactions_collection.update_one(dir_path, {"$set": {'qmmm': 'Already insert'}}, True)

"""
Check barrier which do not have reactant energy
"""

def select_barrier_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'barrier': 'need check'}
    targets = list(qm_collection.find(query))
    return targets

def check_barrier(qm_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_barrier_target(qm_collection)
    # 2. check the job pbs status
    for target in targets:
        if target['irc_equal'] == 'forward equal to reactant':
            reactant_energy = float(target['irc_forward_opt_energy'])
        elif target['irc_equal'] == 'backward equal to reactant':
            reactant_energy = float(target['irc_backward_opt_energy'])
        ts_energy = float(target['ts_energy'])
        barrier = (ts_energy - reactant_energy) * 627.5095
        update_field = {'barrier': barrier}
        qm_collection.update_one(target, {"$set": update_field}, True)

"""
QMMM check.
"""

def select_qmmm_opt_finished_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$and':
             [
                 {"qmmm_opt_reactant_status":
                  {"$in":
                   ['job_success']}},
                 {'qmmm_opt_product_status':
                  {'$in':
                   ['job_success']}}
             ]}
    targets = list(qm_collection.find(query))
    return targets

def check_qmmm_opt_content(dir_path:str, direction:str='reactant') -> str:
    qmmm_reactant_dir = path.join(dir_path, 'QMMM_REACTANT/')
    qmmm_product_dir = path.join(dir_path, 'QMMM_PRODUCT/')
    if direction == 'reactant':
        output_path = path.join(qmmm_reactant_dir, 'qmmm_opt.out')
    else:
        output_path = path.join(qmmm_product_dir, 'qmmm_opt.out')
    try:
        q = QChem(outputfile=output_path)
        if direction == 'reactant':
            reactant = path.join(qmmm_reactant_dir, 'qmmm_opt.xyz')
            q.create_qmmm_geomtry(file_path = reactant)
        else:
            product = path.join(qmmm_product_dir, 'qmmm_opt.xyz')
            q.create_qmmm_geomtry(file_path = product)
        return 'job_success'
    except:
        return 'job_fail'

def check_qmmm_opt_jobs(qm_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_opt_reactant')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_opt_reactant_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_qmmm_opt_content(target['path'], direction='reactant')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_opt_reactant_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_opt_reactant_status': new_status
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_opt_reactant_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_opt_reactant_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='qmmm_opt_product')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_opt_product_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_qmmm_opt_content(target['path'], direction='product')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_opt_product_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_opt_product_status': new_status
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_opt_product_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_opt_product_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    finished_targets = select_qmmm_opt_finished_target(qm_collection)
    for target in finished_targets:
        update_field = {
            'qmmm_opt_status': 'job_success', 'qmmm_freq_opt_status':'job_unrun'
            }
        qm_collection.update_one(target, {"$unset": {'qmmm_opt_reactant_status': '', 'qmmm_opt_product_status': '',
                                                     'qmmm_opt_reactant_jobid': '', 'qmmm_opt_product_jobid': ''}, "$set": update_field}, True)

"""
QMMM FREQ OPT
"""

def select_qmmm_freq_opt_finished_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$and':
             [
                 {"qmmm_freq_opt_reactant_status":
                  {"$in":
                   ['job_success']}},
                 {'qmmm_freq_opt_product_status':
                  {'$in':
                   ['job_success']}}
             ]}
    targets = list(qm_collection.find(query))
    return targets

def check_qmmm_freq_opt_content(dir_path:str, direction:str='reactant', restart:bool=False) -> str:
    qmmm_reactant_dir = path.join(dir_path, 'QMMM_REACTANT/')
    qmmm_product_dir = path.join(dir_path, 'QMMM_PRODUCT/')
    if direction == 'reactant':
        output_path = path.join(qmmm_reactant_dir, 'qmmm_freq_opt.out')
    else:
        output_path = path.join(qmmm_product_dir, 'qmmm_freq_opt.out')
    try:
        q = QChem(outputfile=output_path)
        if direction == 'reactant':
            reactant = path.join(qmmm_reactant_dir, 'qmmm_freq_opt.xyz')
            q.create_qmmm_geomtry(file_path = reactant)
            if restart:
                return 'job_success'
            else:
                return 'restart'
        else:
            product = path.join(qmmm_product_dir, 'qmmm_freq_opt.xyz')
            q.create_qmmm_geomtry(file_path = product)
            if restart:
                return 'job_success'
            else:
                return 'restart'
    except:
        return 'job_fail'

def check_qmmm_freq_opt_jobs(qm_collection:object, restart_times:int = 2):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_freq_opt_reactant')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_freq_opt_reactant_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            times = target['qmmm_freq_opt_reactant_restart_times']
            if times >= restart_times:
                new_status = check_qmmm_freq_opt_content(target['path'], direction='reactant', restart=True)
            else:
                new_status = check_qmmm_freq_opt_content(target['path'], direction='reactant')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_freq_opt_reactant_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_freq_opt_reactant_status': new_status
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_freq_opt_reactant_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_freq_opt_reactant_status': new_status, 'qmmm_freq_opt_reactant_restart_times':times + 1
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='qmmm_freq_opt_product')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_freq_opt_product_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            times = target['qmmm_freq_opt_product_restart_times']
            if times >= restart_times:
                new_status = check_qmmm_freq_opt_content(target['path'], direction='product', restart=True)
            else:
                new_status = check_qmmm_freq_opt_content(target['path'], direction='product')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_freq_opt_product_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_freq_opt_product_status': new_status
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_freq_opt_product_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_freq_opt_product_status': new_status, 'qmmm_freq_opt_product_restart_times':times + 1
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    finished_targets = select_qmmm_freq_opt_finished_target(qm_collection)
    for target in finished_targets:
        update_field = {
            'qmmm_freq_opt_status': 'job_success', 'qmmm_freq_status': 'job_unrun'
            }
        qm_collection.update_one(target, {"$unset": {'qmmm_freq_opt_reactant_status': '', 'qmmm_freq_opt_product_status': '',
                                                     'qmmm_freq_opt_reactant_jobid': '', 'qmmm_freq_opt_product_jobid': ''}, "$set": update_field}, True)

"""
QMMM FREQ
"""

def select_qmmm_freq_finished_target() -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    query = {'$and':
             [
                 {"qmmm_freq_reactant_status":
                  {"$in":
                   ['job_success']}},
                 {'qmmm_freq_product_status':
                  {'$in':
                   ['job_success']}}
             ]
             }
    targets = list(qm_collection.find(query))
    return targets

def check_qmmm_freq_content(dir_path:str, direction:str='reactant') -> Union[str, float]:
    qmmm_reactant_dir = path.join(dir_path, 'QMMM_REACTANT/')
    qmmm_product_dir = path.join(dir_path, 'QMMM_PRODUCT/')

    if direction == 'reactant':
        output_path = path.join(qmmm_reactant_dir, 'qmmm_freq_opt.out')
        reactant_path = path.join(qmmm_reactant_dir, 'qmmm_final.xyz')
    else:
        output_path = path.join(qmmm_product_dir, 'qmmm_freq_opt.out')
        reactant_path = path.join(qmmm_product_dir, 'qmmm_final.xyz')

    try:
        q = QChem(outputfile=output_path)
        freqs = q.get_frequencies()
        nnegfreq = sum(1 for freq in freqs if freq < 0.0)
        if nnegfreq > 0:
            return 'Have negative frequency', 0.0
        else:
            q.create_qmmm_geomtry(reactant_path)
            energy = q.get_energy()
            zpe = q.get_zpe()
            energy += zpe
            return 'job_success', energy
    except:
        return 'job_fail', 0.0

def check_qmmm_freq_jobs(qm_collection:object, reactions_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_freq_reactant')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_freq_reactant_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_freq_content(target['path'], direction='reactant')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_freq_reactant_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_freq_reactant_status': new_status, 'qmmm_freq_reactant_energy':energy
                    }
                update_field_reaction = {
                    'qmmm_freq_reactant_status': new_status, 'qmmm_freq_reactant_energy':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_freq_reactant_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_freq_reactant_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_freq_reactant_status': new_status
                    }
                update_field_reaction = {
                        'qmmm_freq_reactant_status': new_status, 'qmmm_freq_reactant_energy':energy
                    }
            reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
            reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='qmmm_freq_product')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_freq_product_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_freq_content(target['path'], direction='product')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_freq_product_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_freq_product_status': new_status, 'qmmm_freq_product_energy':energy
                    }
                update_field_reaction = {
                    'qmmm_freq_product_status': new_status, 'qmmm_freq_product_energy':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_freq_product_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_freq_product_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_freq_product_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_freq_product_status': new_status, 'qmmm_freq_product_energy':energy
                    }
            reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
            reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
            qm_collection.update_one(target, {"$set": update_field}, True)

    finished_targets = select_qmmm_freq_finished_target()
    for target in finished_targets:
        reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
        update_field = {
            'qmmm_freq_status': 'job_success', 'qmmm_refine_status':'job_unrun'
            }
        # qmmm_delta_H = target['qmmm_freq_product_energy'] - target['qmmm_freq_reactant_energy']
        # update_field_reaction = {'qmmm_delta_H':qmmm_delta_H}
        # reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
        qm_collection.update_one(target, {"$unset": {'qmmm_freq_reactant_status': '', 'qmmm_freq_product_status': '',
                                                     'qmmm_freq_reactant_jobid': '', 'qmmm_freq_product_jobid': ''}, "$set": update_field}, True)

"""
QMMM FREQ TS
"""

def check_qmmm_freq_ts_content(dir_path:str, restart:bool=False) -> str:
    qmmm_ts_dir = path.join(dir_path, 'QMMM_TS/')
    output_path = path.join(qmmm_ts_dir, 'qmmm_freq_ts.out')
    ts = path.join(qmmm_ts_dir, 'qmmm_ts.xyz')
    try:
        q = QChem(outputfile=output_path)
        q.create_qmmm_geomtry(file_path = ts)
        if restart:
            return 'job_success'
        else:
            return 'restart'
    except:
        return 'job_fail'

def check_qmmm_freq_ts_jobs(qm_collection:object, restart_times:int = 2):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_freq_ts')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_freq_ts_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            times = target['qmmm_freq_ts_restart_times']
            if times >= restart_times:
                new_status = check_qmmm_freq_ts_content(target['path'], restart=True)
            else:
                new_status = check_qmmm_freq_ts_content(target['path'])
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_freq_ts_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_freq_ts_status': new_status, "qmmm_ts_freq_status": "job_unrun"
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_freq_ts_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_freq_ts_status': new_status, 'qmmm_freq_ts_restart_times': times + 1
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
QMMM TS FREQ
"""

def check_qmmm_ts_freq_content(dir_path:str) -> Union[str, float]:
    qmmm_ts_dir = path.join(dir_path, 'QMMM_TS/')
    output_path = path.join(qmmm_ts_dir, 'qmmm_freq.out')
    ts_path = path.join(qmmm_ts_dir, 'qmmm_final.xyz')

    try:
        q = QChem(outputfile=output_path)
        freqs = q.get_frequencies()
        nnegfreq = sum(1 for freq in freqs if freq < 0.0)
        if nnegfreq > 1:
            return "Have more than one imaginary frequency", 0.0
        elif nnegfreq == 0:
            return "All positive frequency", 0.0
        else:
            q.create_qmmm_geomtry(ts_path)
            energy = q.get_energy()
            zpe = q.get_zpe()
            energy += zpe
            return 'job_success', energy
    except:
        return 'job_fail', 0.0

def check_qmmm_ts_freq_jobs(qm_collection:object, reactions_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_ts_freq')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_ts_freq_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_ts_freq_content(target['path'])
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_ts_freq_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_ts_freq_status': new_status, "qmmm_ts_energy": energy, 'qmmm_ts_refine': 'job_unrun'
                    }
                update_field_reaction = {
                    'qmmm_ts_freq_status': new_status, 'qmmm_ts_energy':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_ts_freq_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_ts_freq_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_ts_freq_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_ts_freq_status': new_status, 'qmmm_ts_energy':energy
                    }
            reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
            reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
            qm_collection.update_one(target, {"$set": update_field}, True)

"""
QMMM REFINE
"""

def select_qmmm_refine_finished_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$and':
             [
                 {"qmmm_refine_reactant_status":
                  {"$in":
                   ['job_success']}},
                 {'qmmm_refine_product_status':
                  {'$in':
                   ['job_success']}},
                 {"qmmm_refine_ts_status":
                  {"$in":
                   ['job_success']}}
             ]}
    targets = list(qm_collection.find(query))
    return targets

def check_qmmm_refine_content(dir_path:str, direction:str='reactant') -> Union[str, float]:
    qmmm_reactant_dir = path.join(dir_path, 'QMMM_REACTANT/')
    qmmm_product_dir = path.join(dir_path, 'QMMM_PRODUCT/')
    qmmm_ts_dir = path.join(dir_path, 'QMMM_TS/')
    if direction == 'reactant':
        output_path = path.join(qmmm_reactant_dir, 'qmmm_sp.out')
    elif direction == 'product':
        output_path = path.join(qmmm_product_dir, 'qmmm_sp.out')
    else:
        output_path = path.join(qmmm_ts_dir, 'qmmm_sp.out')

    try:
        q = QChem(outputfile=output_path)
        energy = q.get_energy()
        return 'job_success', energy
    except:
        return 'job_fail', 0.0

def check_qmmm_refine_jobs(qm_collection:object, reactions_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_targets(qm_collection, job_name='qmmm_refine_reactant')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_refine_reactant_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_refine_content(target['path'], direction='reactant')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_refine_reactant_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_refine_reactant_status': new_status, 'qmmm_sp_reactant':energy
                    }
                update_field_reaction = {
                    'qmmm_refine_reactant_status': new_status, 'qmmm_sp_reactant':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_refine_reactant_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_refine_reactant_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_refine_reactant_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_refine_reactant_status': new_status
                    }
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='qmmm_refine_product')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_refine_product_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_refine_content(target['path'], direction='product')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_refine_product_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_refine_product_status': new_status, 'qmmm_sp_product':energy
                    }
                update_field_reaction = {
                    'qmmm_refine_product_status': new_status, 'qmmm_sp_product':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_refine_product_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_refine_product_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_refine_product_status': new_status
                    }
                update_field_reaction = {
                    'qmmm_refine_product_status': new_status
                    }
            reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
            reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
            qm_collection.update_one(target, {"$set": update_field}, True)

    targets = select_targets(qm_collection, job_name='qmmm_refine_ts')
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['qmmm_refine_ts_jobid']
        new_status = check_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_qmmm_refine_content(target['path'], direction='ts')
        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['qmmm_refine_ts_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                    'qmmm_refine_ts_status': new_status, 'qmmm_sp_ts':energy
                    }
            elif new_status == "job_running" or new_status == "job_queueing" or new_status == "job_launched":
                update_field = {
                    'qmmm_refine_ts_status': new_status
                    }
            else:
                update_field = {
                    'qmmm_refine_ts_status': new_status
                    }
            reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
            reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
            qm_collection.update_one(target, {"$set": update_field}, True)

    finished_targets = select_qmmm_refine_finished_target(qm_collection)
    for target in finished_targets:
        delta_H = (target['qmmm_sp_product'] - target['qmmm_sp_reactant']) * 627.5095
        barrier = (target['qmmm_sp_ts'] - target['qmmm_sp_reactant']) * 627.5095
        update_field = {
            'qmmm_refine_status': 'job_success', 'qmmm_delta_H':delta_H, 'qmmm_barrier':barrier
            }
        update_field_reaction = {
            'qmmm_refine_status': 'job_success', 'qmmm_delta_H':delta_H, 'qmmm_barrier':barrier
            }
        reaction_target = list(reactions_collection.find({'path':target['path']}))[0]
        reactions_collection.update_one(reaction_target, {"$set": update_field_reaction}, True)
        qm_collection.update_one(target, {"$unset": {'qmmm_refine_reactant_status': '', 'qmmm_refine_reactant_jobid': '',
                                                     'qmmm_refine_product_status': '', 'qmmm_refine_product_jobid': '',
                                                     'qmmm_refine_ts_status': '', 'qmmm_refine_ts_jobid': ''}, 
                                                     "$set": update_field}, True)

"""
QMMM side fail check.
If qmmm ts fail then delete the reactant and product opt (while the other side still running).
"""

def select_qmmm_freq_ts_side_fail_target(qm_collection:object) -> list:
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    query = {'$or':[{'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_freq_opt_reactant_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            },
            {'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_freq_opt_product_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            },
            {'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_opt_reactant_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            },
            {'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_opt_product_status':
                  {'$in':
                    ["job_running", "job_queueing", "job_launched"]}}
            ]
            },
            {'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_freq_opt_status':
                  {'$in':
                    ["job_unrun"]}}
            ]
            },
            {'$and':
            [
                {"qmmm_freq_ts_status":
                {"$in":
                    ['job_fail']}},
                {'qmmm_opt_status':
                  {'$in':
                    ["job_unrun"]}}
            ]
            }]}

    targets = list(qm_collection.find(query))
    return targets

def check_qmmm_freq_ts_side_fail_jobs(qm_collection:object):
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_qmmm_freq_ts_side_fail_target(qm_collection)
    # 2. check the job pbs_status
    for target in targets:
        if target['qmmm_freq_opt_reactant_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['qmmm_freq_opt_reactant_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'qmmm_freq_opt_reactant_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['qmmm_freq_opt_product_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['qmmm_freq_opt_product_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'qmmm_freq_opt_reactant_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['qmmm_freq_opt_status'] in ["job_unrun"]:
            update_field = {
                'qmmm_freq_opt_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['qmmm_opt_reactant_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['qmmm_opt_reacrant_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'qmmm_opt_reactant_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['qmmm_opt_product_status'] in ["job_running", "job_queueing", "job_launched"]:
            job_id = target['qmmm_opt_product_jobid']
            commands = ['qdel', job_id]
            process = subprocess.Popen(commands,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            update_field = {
                'qmmm_opt_product_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['qmmm_opt_status'] in ["job_unrun"]:
            update_field = {
                'qmmm_opt_status': 'qmmm ts fail'
                }
            qm_collection.update_one(target, {"$set": update_field}, True)

def check_jobs(refine=True, cluster_bond_path=None, level_of_theory='ORCA'):
    qm_collection = db['qm_calculate_center']
    reactions_collection = db['reactions']
    statistics_collection = db['statistics']
    config_collection = db['config']
    if cluster_bond_path:
        # use the checker.py path as the reference
        checker_path = os.path.realpath(sys.argv[0])
        ard_path = os.path.dirname(os.path.dirname(checker_path))
        cluster_bond_path = path.join(ard_path, 'script/bonds.txt')
        fixed_atom_path = path.join(ard_path, 'script/fixed_atom.txt')
    
    # If the ssm perform by orca with xtb GFN2-xtb, then refine the TS is a good choice.  Get a better initial guess
    check_ssm_jobs(qm_collection, refine=refine, thershold = 80.0)  # TS guess energy filter
    check_ts_refine_jobs(qm_collection, threshold = -50.0)  # Imaginary freq should smaller than threshold
    check_ts_jobs(qm_collection, threshold = -50.0) # Imaginary freq should smaller than threshold
    check_irc_jobs(qm_collection)
    check_irc_opt_jobs(qm_collection, level_of_theory=level_of_theory)
    check_irc_opt_side_fail_jobs(qm_collection)
    check_irc_equal(qm_collection, cluster_bond_path = cluster_bond_path, fixed_atom_path = fixed_atom_path)
    check_barrier(qm_collection)
    insert_reaction(qm_collection, reactions_collection)
    insert_ard(qm_collection, reactions_collection, statistics_collection, config_collection, barrier_threshold=60.0)

    check_qmmm_opt_jobs(qm_collection)
    check_qmmm_freq_opt_jobs(qm_collection, restart_times = 2)
    check_qmmm_freq_ts_jobs(qm_collection, restart_times = 2)
    check_qmmm_freq_jobs(qm_collection, reactions_collection)
    check_qmmm_ts_freq_jobs(qm_collection, reactions_collection)
    check_qmmm_refine_jobs(qm_collection, reactions_collection)
    check_qmmm_freq_ts_side_fail_jobs(qm_collection)

check_jobs(refine=True, cluster_bond_path=True, level_of_theory='ORCA')
