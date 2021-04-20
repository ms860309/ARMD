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


def select_ard_target():
    # when lower generation finished return target else return []
    qm_collection = db['qm_calculate_center']
    reg_query = {'ard_status':
                 {"$in":
                  ["job_unrun"]
                  }
                 }
    targets = list(qm_collection.find(reg_query))
    return targets


def launch_ard_jobs(ncpus=8, mpiprocs=1, ompthreads=8):

    qm_collection = db['qm_calculate_center']
    pool_collection = db['pool']
    status_collection = db['status']

    if pool_collection.estimated_document_count() == 0:
        print(highlight_text('Starting ARD network exploring'))
        script_path = path.join(path.dirname(
            path.dirname(path.abspath(__file__))), 'script')
        if os.path.exists(script_path):
            os.chdir(script_path)
        subfile = create_ard_sub_file(
            script_path, script_path, 1, 'reactant.xyz', ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        # first reactant need to add to pool
        initial_reactant = next(pybel.readfile(
            'xyz', path.join(script_path, 'reactant.xyz')))
        initial_reactant_inchi_key = initial_reactant.write('inchiKey').strip()
        pool_collection.insert_one(
            {'reactant_inchi_key': initial_reactant_inchi_key, 'generations': 1})
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
        status_collection.insert_one(
            {'status': 'ARD had launched', 'jobid': job_id})
    else:
        targets = select_ard_target()
        use_irc_query = {'reactant_smiles': 'initial reactant'}
        use_irc = list(qm_collection.find(use_irc_query))[0]['use_irc']
        count = 0
        for target in targets:
            if use_irc == '0':
                dir_path, gen_num, ard_ssm_equal = target['path'], target['generations'], target['ard_ssm_equal']
                script_path = path.join(path.dirname(
                    path.dirname(dir_path)), 'script')
                if ard_ssm_equal == 'not_equal':
                    next_reactant = 'ssm_product.xyz'
                    subfile = create_ard_sub_file(
                        dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
                else:
                    next_reactant = 'product.xyz'
                    subfile = create_ard_sub_file(
                        dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
            else:
                dir_path, gen_num = target['path'], target['generations']
                next_reactant = path.join(dir_path, 'irc_reactant.xyz')
                script_path = path.join(path.dirname(
                    path.dirname(dir_path)), 'script')
                subfile = create_ard_sub_file(
                    dir_path, script_path, gen_num + 1, next_reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)

            os.chdir(dir_path)

            cmd = 'qsub {}'.format(subfile)
            process = subprocess.Popen([cmd],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, shell=True)
            stdout, stderr = process.communicate()
            # get job id from stdout, e.g., "106849.h81"
            job_id = stdout.decode().replace("\n", "")
            # update status job_launched
            update_ard_status(target, job_id)
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
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(
        ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(script_path)
    nes1 = 'source ~/.bashrc\nexport MKL_NUM_THREADS={}\nexport OMP_NUM_THREADS={}\nexport OMP_STACKSIZE={}G\n'.format(
        ncpus, ompthreads, mem)
    nes2 = 'conda activate ard'
    command = 'python {} {} {} -bonds {} -constraint {} -fixed_atom {} -generations {}'.format(
        ard_path, input_path, product_xyz_path, bonds_path, constraint, fixed_atom, gen_num)
    deactivate = 'conda deactivate'

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(
            shell, pbs_setting, target_path, nes1, nes2, command, deactivate))

    return subfile


def update_ard_status(target, job_id):
    qm_collection = db['qm_calculate_center']
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



launch_ard_jobs(ncpus=1, mpiprocs=1, ompthreads=1)
