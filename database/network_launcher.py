from openbabel import pybel
from connect import db
import subprocess
import os
from os import path
import shutil
import sys
from helper import highlight_text, print_header
sys.path.append(path.join(path.dirname(
    path.dirname(path.abspath(__file__))), 'script'))

# Manual run 0th generations


def select_ard_target(qm_collection:object):
    # when lower generation finished return target else return []
    reg_query = {'ard_status':
                 {"$in":
                  ["job_unrun"]
                  }
                 }
    targets = list(qm_collection.find(reg_query))
    return targets


def launch_ard_jobs(qm_collection:object, config_collection:object, status_collection:object, ncpus:int=1, mpiprocs:int=1, ompthreads:int=1, qmmm:bool=True):

    if config_collection.estimated_document_count() == 0:
        print(highlight_text('Starting ARD network exploring'))
        script_path = path.join(path.dirname(path.dirname(path.abspath(__file__))), 'script')
        if path.exists(script_path):
            os.chdir(script_path)
        subfile = create_ard_sub_file(script_path, script_path, 1, 'reactant.xyz', ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        # first reactant need to add to config
        initial_reactant = next(pybel.readfile('xyz', path.join(script_path, 'reactant.xyz')))
        initial_reactant_inchi_key = initial_reactant.write('inchiKey').strip()
        config_collection.insert_one({'reactant_inchi_key': initial_reactant_inchi_key, 'generations': 1})
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        print('ARD had launched')
        print('jobid is {}'.format(job_id))
        status_collection.insert_one({'status': 'ARD had launched', 'jobid': job_id})
    else:
        targets = select_ard_target(qm_collection)
        use_irc = list(config_collection.find({'generations':1}))[0]['use_irc']
        count = 0
        for target in targets:
            if use_irc == '0':
                dir_path, gen_num, ard_ssm_equal = target['path'], target['generations'], target['ard_ssm_equal']
                script_path = path.join(path.dirname(
                    path.dirname(dir_path)), 'script')
                if ard_ssm_equal == 'not_equal':
                    next_reactant = 'ssm_product.xyz'
                    subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
                else:
                    next_reactant = 'product.xyz'
                    subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
            else:
                dir_path, gen_num = target['path'], target['generations']
                qmmm_product_path = path.join(dir_path, 'QMMM_PRODUCT/qmmm_final.xyz')
                next_reactant_path = path.join(dir_path, 'qmmm_reactant.xyz')
                if qmmm:
                    shutil.copyfile(qmmm_product_path, next_reactant_path)
                    next_reactant = 'qmmm_reactant.xyz'
                else:
                    next_reactant = 'irc_reactant.xyz'
                script_path = path.join(path.dirname(path.dirname(dir_path)), 'script')
                subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)

            os.chdir(dir_path)

            cmd = 'qsub {}'.format(subfile)
            process = subprocess.Popen([cmd],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, shell=True)
            stdout, stderr = process.communicate()
            # get job id from stdout, e.g., "106849.h81"
            job_id = stdout.decode().replace("\n", "")
            # update status job_launched
            update_ard_status(qm_collection, target, job_id)
            count += 1
        print(highlight_text('ARD'))
        print('\nARD launced {} jobs\n'.format(count))

def create_ard_sub_file(dir_path, script_path, gen_num, next_reactant, ncpus=4, mpiprocs=1, ompthreads=4, mem=1):
    subfile = path.join(dir_path, 'ard.job')
    product_xyz_path = path.join(dir_path, next_reactant)
    ard_path = path.join(script_path, 'ard.py')
    input_path = path.join(script_path, 'input.txt')
    bonds_path = path.join(script_path, 'bonds.txt')
    constraint = path.join(script_path, 'constraint.txt')
    fixed_atom = path.join(script_path, 'fixed_atom.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(script_path)
    nes1 = 'source ~/.bashrc\nexport MKL_NUM_THREADS={}\nexport OMP_NUM_THREADS={}\nexport OMP_STACKSIZE={}G\n'.format(ncpus, ompthreads, mem)
    nes2 = 'conda activate ard'
    command = 'python {} {} {} -bonds {} -constraint {} -fixed_atom {} -generations {}'.format(ard_path, input_path, product_xyz_path, bonds_path, constraint, fixed_atom, gen_num)
    deactivate = 'conda deactivate'

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, command, deactivate))

    return subfile

def update_ard_status(qm_collection, target, job_id):
    update_field = {"ard_status": "job_launched", "ard_jobid": job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)

def check_ard_job_status(job_id):
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """

    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if 'Unknown Job Id' in stderr.decode():
        return 'off_queue'
    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return 'job_running'
    elif stdout[idx+2] == 'Q':
        return 'job_queueing'
    else:
        return 'job_launched'



qm_collection = db['qm_calculate_center']
config_collection = db['config']
status_collection = db['status']

target = list(config_collection.find({'generations': 1}))[0]
qmmm = target['use_qmmm']
if qmmm == '1':
    qmmm = True
else:
    qmmm = False
launch_ard_jobs(qm_collection, config_collection, status_collection, ncpus=1, mpiprocs=1, ompthreads=1, qmmm=qmmm)
