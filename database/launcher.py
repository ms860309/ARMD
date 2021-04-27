from connect import db
from helper import highlight_text, print_header
import subprocess
import os
from os import path
import shutil
from shutil import copyfile
import sys
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'script'))

class LaunchError(Exception):
    pass

"""
Global utils
"""

def select_targets(qm_collection:object, job_name:str) -> list:
    keyname = '{}_status'.format(job_name)
    query = {keyname:
             {"$in":
              ["job_unrun", "restart"]
              }}
    targets = list(qm_collection.find(query))
    return targets


def update_status(qm_collection:object, target:object, job_name:str, job_id:str):
    keyname = '{}_status'.format(job_name)
    keyjobid = '{}_jobid'.format(job_name)
    update_field = {keyname: "job_launched", keyjobid: job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)

"""
Submmit SSM calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def launch_ssm_jobs(qm_collection:object, num:int=100, level_of_theory:str='QCHEM', ncpus:int=1, mpiprocs:int=1, ompthreads:int=1):
    targets = select_targets(qm_collection, job_name='ssm')
    count = 0
    for target in targets[:num]:
        SSM_dir_path = path.join(target['path'], 'SSM/')
        os.mkdir(SSM_dir_path)
        os.chdir(SSM_dir_path)
        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_ssm_sub_file(target['path'], SSM_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        elif level_of_theory.upper() == 'ORCA':
            subfile = create_orca_ssm_sub_file(target['path'], SSM_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_status(qm_collection, target, job_name = 'ssm', job_id = job_id)
        count += 1
    print(highlight_text('SSM'))
    print('\nSSM launced {} jobs\n'.format(count))

def create_qchem_ssm_sub_file(dir_path:str, SSM_dir_path:str, ncpus:int=1, mpiprocs:int=1, ompthreads:int=1) -> str:
    subfile = path.join(SSM_dir_path, 'ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'qchem_qstart')
    ssm_args = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'ssm_argument')
    frozen_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'frozen.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    nes1 = 'module load qchem'
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nconda activate ard'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    if os.path.exists(frozen_file):
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} -frozen_coord_idx_file {} '.format(xyz_file, isomers, lot_inp_file, frozen_file)
    else:
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} '.format(xyz_file, isomers, lot_inp_file)
    with open(ssm_args, 'r') as f:
        lines = f.read().splitlines()
    command = command + ' '.join(lines) + ' > status.log 2>&1 '
    clean_scratch = """rm -r $QCSCRATCH\nconda deactivate\necho "It took $(($(date +'%s') - $start)) seconds" """
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, scratch, command, clean_scratch))

    return subfile

def create_orca_ssm_sub_file(dir_path:str, SSM_dir_path:str, ncpus:int=1, mpiprocs:int=1, ompthreads:int=1, mem:int=1) -> str:
    subfile = path.join(SSM_dir_path, 'ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    ssm_args = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'ssm_argument')
    frozen_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'frozen.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nexport MKL_NUM_THREADS={}\nexport OMP_NUM_THREADS={}\nexport OMP_STACKSIZE={}G\nconda activate ard\n'.format(ncpus, ompthreads, mem)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    if os.path.exists(frozen_file):
        command = 'gsm -xyzfile {} -mode SE_GSM -package xTB_lot -isomers {} -frozen_coord_idx_file {} '.format(xyz_file, isomers, frozen_file)
    else:
        command = 'gsm -xyzfile {} -mode SE_GSM -package xTB_lot -isomers {} '.format(xyz_file, isomers)

    with open(ssm_args, 'r') as f:
        lines = f.read().splitlines()
    command = command + ' '.join(lines) + ' > status.log 2>&1 '
    clean_scratch = """rm -r $QCSCRATCH\nconda deactivate\necho "It took $(($(date +'%s') - $start)) seconds" """
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes2, scratch, command, clean_scratch))

    return subfile

"""
Submmit TS refine calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
Here refine is not to get a high level TS.
It's just for ORCA with xtb GFN2-xtb level of theory SSM. To get a better TS initial guess
"""

def launch_ts_refine_jobs(qm_collection:object, num:int=100, ncpus:int=1, mpiprocs:int=1, ompthreads:int=1):
    targets = select_targets(qm_collection, job_name='ts_refine')
    count = 0
    for target in targets[:num]:
        TS_dir_path = path.join(target['path'], 'TS/')
        os.mkdir(TS_dir_path)
        os.chdir(TS_dir_path)

        SSM_dir_path = path.join(target['path'], 'SSM/')
        subfile = create_ts_refine_sub_file(SSM_dir_path, TS_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_status(qm_collection, target, job_name = 'ts_refine', job_id = job_id)
        count += 1
    print(highlight_text('TS refine'))
    print('\nTS refine launced {} jobs\n'.format(count))

def create_ts_refine_sub_file(SSM_dir_path:str, TS_dir_path:str, ncpus:int=1, mpiprocs:int=1, ompthreads:int=1, mem:int=2) -> str:
    tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts_refine.in')
    subfile = path.join(TS_dir_path, 'ts_refine.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'orca_refine_ts')
    with open(ts_lot) as f:
        config = [line.strip() for line in f]
    path2orca = os.popen('which orca').read().rstrip()
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\nmodule load orca\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export MKL_NUM_THREADS={}\nexport OMP_NUM_THREADS={}\nexport OMP_STACKSIZE={}G\nexport QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/ts_refine.in $QCSCRATCH/\ncd $QCSCRATCH\n'.format(ncpus, ompthreads, mem)
    command = '{} $QCSCRATCH/ts_refine.in >> $PBS_O_WORKDIR/ts_refine.out'.format(path2orca)
    copy_the_refine_xyz = 'cp $QCSCRATCH/ts_refine.xyz $PBS_O_WORKDIR'
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('\n%pal\nnprocs {}\nend\n\n'.format(ncpus))
        f2.write('*xyz 0 1\n')

        for line in lines[2:]:
            f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, copy_the_refine_xyz, clean_scratch))

    return subfile

"""
Submmit TS calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def launch_ts_jobs(qm_collection:object, num:int=100, level_of_theory:str='ORCA', ncpus:int=4, mpiprocs:int=4, ompthreads:int=1, Hcap:int=None):
    targets = select_targets(qm_collection, job_name='ts')
    count = 0
    for target in targets[:num]:
        SSM_dir_path = path.join(target['path'], 'SSM/')
        TS_dir_path = path.join(target['path'], 'TS/')
        if os.path.exists(TS_dir_path):
            os.chdir(TS_dir_path)
        else:
            os.mkdir(TS_dir_path)
            os.chdir(TS_dir_path)

        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        elif level_of_theory == 'ORCA':
            subfile = create_orca_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads, Hcap=Hcap)
        else:
            raise LaunchError('Unsupported level of theory')
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_status(qm_collection, target, job_name = 'ts', job_id = job_id)
        count += 1
    print(highlight_text('TS'))
    print('\nTS launced {} jobs\n'.format(count))

def create_qchem_ts_sub_file(SSM_dir_path:str, TS_dir_path:str, ncpus:int=4, mpiprocs:int=1, ompthreads:int=4) -> str:
    refine_path = path.join(TS_dir_path, 'ts_refine.xyz')
    if os.path.exists(refine_path):
        tsnode_path = refine_path
    else:
        tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts.in')
    ts_output_file = path.join(TS_dir_path, 'ts.out')
    subfile = path.join(TS_dir_path, 'ts.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'qchem_freq_ts_freq.lot')
    with open(ts_lot) as f:
        config = [line.strip() for line in f]
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(TS_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, ts_input_file, ts_output_file)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """
    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))

    return subfile

def create_orca_ts_sub_file(SSM_dir_path:str, TS_dir_path:str, ncpus:int=4, mpiprocs:int=4, ompthreads:int=1, Hcap:int=None) -> str:
    refine_path = path.join(TS_dir_path, 'ts_refine.xyz')
    if os.path.exists(refine_path):
        tsnode_path = refine_path
    else:
        tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts_geo.in')
    subfile = path.join(TS_dir_path, 'ts.job')
    base_dir_path = path.join(path.dirname(path.dirname(
        path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'orca_freq_ts_freq.lot')
    with open(ts_lot) as f:
        config = [line.strip() for line in f]
    path2orca = os.popen('which orca').read().rstrip()
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\nmodule load orca\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/ts_geo.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command = '{} $QCSCRATCH/ts_geo.in >> $PBS_O_WORKDIR/ts_geo.out'.format(path2orca)
    copy_the_refine_xyz = 'cp $QCSCRATCH/ts_geo.xyz $PBS_O_WORKDIR'
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('\n%pal\nnprocs {}\nend\n\n'.format(ncpus))
        f2.write('*xyz 0 1\n')
        if Hcap:
            num = len(lines[2:]) - Hcap
            for idx, line in enumerate(lines[2:]):
                if idx < num:
                    f2.write(line + '\n')
                else:
                    f2.write(line + ' m=100000000000000000' + '\n')
        else:
            for line in lines[2:]:
                f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, copy_the_refine_xyz, clean_scratch))

    return subfile

"""
Submmit IRC(Intrinsic Reaction Coordinate, in Qchem called 'rpath') calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def launch_irc_jobs(qm_collection:object, num:int=100, ncpus:int=8, mpiprocs:int=8, ompthreads:int=1):
    targets = select_targets(qm_collection, job_name='irc')
    count = 0
    for target in targets[:num]:
        IRC_dir_path = path.join(target['path'], 'IRC/')
        os.mkdir(IRC_dir_path)
        os.chdir(IRC_dir_path)

        TS_dir_path = path.join(target['path'], 'TS/')
        subfile = create_irc_sub_file(
            TS_dir_path, IRC_dir_path, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_status(qm_collection, target, job_name = 'irc', job_id = job_id)
        count += 1
    print(highlight_text('IRC'))
    print('\nIRC launced {} jobs\n'.format(count))

def create_irc_sub_file(TS_dir_path:str, IRC_dir_path:str, ncpus:int=8, mpiprocs:int=8, ompthreads:int=1) -> str:
    ts_geo_path = path.join(TS_dir_path, 'ts_geo.xyz')
    irc_input_file = path.join(IRC_dir_path, 'pysisyphus_irc.yaml')
    subfile = path.join(IRC_dir_path, 'irc.job')
    new_ts_geo_path = path.join(IRC_dir_path, 'ts_geo.xyz')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(IRC_dir_path)))), 'config')
    irc_lot = path.join(base_dir_path, 'pysisyphus_irc.yaml')
    copyfile(irc_lot, irc_input_file)
    copyfile(ts_geo_path, new_ts_geo_path)
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(IRC_dir_path)
    nes = 'source ~/.bashrc\nmodule load orca\nconda activate ard\nexport TMPDIR=/tmp/$PBS_JOBID\nmkdir -p $TMPDIR\n'
    command = 'pysis pysisyphus_irc.yaml'
    deactivate = 'conda deactivate\nrm -r $TMPDIR'
    clean_scratch = """echo "It took $(($(date +'%s') - $start)) seconds" """
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes, command, deactivate, clean_scratch))

    return subfile

"""
Submmit opt job which is from irc
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def launch_irc_opt_jobs(qm_collection:object, num:int=100, level_of_theory:str='ORCA', ncpus:int=8, mpiprocs:int=8, ompthreads:int=1, Hcap:int=None):
    targets = select_targets(qm_collection, job_name='irc_opt')
    count = 0
    for target in targets[:num]:
        IRC_dir_path = path.join(target['path'], 'IRC/')
        os.chdir(IRC_dir_path)

        first_output = os.path.join(IRC_dir_path, 'finished_first.xyz')
        last_output = os.path.join(IRC_dir_path, 'finished_last.xyz')
        forward_end_output = os.path.join(IRC_dir_path, 'forward_end_opt.xyz')
        backward_end_output = os.path.join(
            IRC_dir_path, 'backward_end_opt.xyz')
        if not os.path.exists(forward_end_output):
            forward_end_output = first_output
        if not os.path.exists(backward_end_output):
            backward_end_output = last_output
        if level_of_theory.upper() == 'QCHEM':
            subfile_1, subfile_2 = create_qchem_irc_opt_sub_file(IRC_dir_path, forward_end_output, backward_end_output, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        elif level_of_theory == 'ORCA':
            subfile_1, subfile_2 = create_orca_irc_opt_sub_file(IRC_dir_path, forward_end_output, backward_end_output, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads, Hcap=Hcap)

        cmd_1 = 'qsub {}'.format(subfile_1)
        process = subprocess.Popen([cmd_1],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_1 = stdout.decode().replace("\n", "")

        cmd_2 = 'qsub {}'.format(subfile_2)
        process = subprocess.Popen([cmd_2],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_2 = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(qm_collection, target, job_id_1, job_id_2)
        count += 1
    print(highlight_text('IRC OPT'))
    print('\nIRC opt launced {} jobs (forward + backward)\n'.format(count * 2))

def create_qchem_irc_opt_sub_file(irc_path:str, forward:str, backward:str, ncpus:int=4, mpiprocs:int=1, ompthreads:int=4) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(irc_path)))), 'config')
    irc_opt_lot = path.join(base_dir_path, 'qchem_opt_freq.lot')
    irc_opt_forward_input = path.join(irc_path, 'irc_forward.in')
    subfile_1 = path.join(irc_path, 'irc_forward_opt.job')
    irc_opt_backward_input = path.join(irc_path, 'irc_backward.in')
    subfile_2 = path.join(irc_path, 'irc_backward_opt.job')

    with open(irc_opt_lot) as f:
        config = [line.strip() for line in f]
    with open(forward, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_forward_input, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')
    with open(backward, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_backward_input, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(irc_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4_1 = 'qchem -nt {} irc_forward.in irc_forward.out'.format(ncpus)
    nes4_2 = 'qchem -nt {} irc_backward.in irc_backward.out'.format(ncpus)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """
    with open(subfile_1, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4_1, clean_scratch))
    with open(subfile_2, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4_2, clean_scratch))

    return subfile_1, subfile_2

def create_orca_irc_opt_sub_file(irc_path:str, forward:str, backward:str, ncpus:int=4, mpiprocs:int=4, ompthreads:int=1, Hcap:int=None) -> None:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(irc_path)))), 'config')
    irc_opt_lot = path.join(base_dir_path, 'orca_opt_freq.lot')
    irc_opt_forward_input = path.join(irc_path, 'irc_forward.in')
    subfile_1 = path.join(irc_path, 'irc_forward_opt.job')
    irc_opt_backward_input = path.join(irc_path, 'irc_backward.in')
    subfile_2 = path.join(irc_path, 'irc_backward_opt.job')
    with open(irc_opt_lot) as f:
        config = [line.strip() for line in f]
    path2orca = os.popen('which orca').read().rstrip()
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\nmodule load orca\n'.format(ncpus, mpiprocs, ompthreads)
    scratch_1 = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/irc_forward.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command_1 = '{} $QCSCRATCH/irc_forward.in >> $PBS_O_WORKDIR/irc_forward.out'.format(path2orca)
    copy_the_refine_xyz_1 = 'cp $QCSCRATCH/irc_forward.xyz $PBS_O_WORKDIR'
    scratch_2 = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/irc_backward.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command_2 = '{} $QCSCRATCH/irc_backward.in >> $PBS_O_WORKDIR/irc_backward.out'.format(path2orca)
    copy_the_refine_xyz_2 = 'cp $QCSCRATCH/irc_backward.xyz $PBS_O_WORKDIR'
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(forward, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_forward_input, 'w') as f2:
        if '$new_job' in config:
            job2idx = config.index('$new_job')
        else:
            job2idx = len(config)
        for line in config[:job2idx]:
            f2.write(line + '\n')
        f2.write('\n%pal\nnprocs {}\nend\n\n'.format(ncpus))
        f2.write('*xyz 0 1\n')
        if Hcap:
            num = len(lines[2:]) - Hcap
            for idx, line in enumerate(lines[2:]):
                if idx < num:
                    f2.write(line + '\n')
                else:
                    f2.write(line + ' m=100000000000000000' + '\n')
        else:
            for line in lines[2:]:
                f2.write(line + '\n')
        f2.write('*\n\n')
        for line in config[job2idx:]:
            f2.write(line + '\n')

    with open(backward, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_backward_input, 'w') as f2:
        if '$new_job' in config:
            job2idx = config.index('$new_job')
        else:
            job2idx = len(config)
        for line in config[:job2idx]:
            f2.write(line + '\n')
        f2.write('\n%pal\nnprocs {}\nend\n\n'.format(ncpus))
        f2.write('*xyz 0 1\n')
        if Hcap:
            num = len(lines[2:]) - Hcap
            for idx, line in enumerate(lines[2:]):
                if idx < num:
                    f2.write(line + '\n')
                else:
                    f2.write(line + ' m=100000000000000000' + '\n')
        else:
            for line in lines[2:]:
                f2.write(line + '\n')
        f2.write('*\n\n')
        for line in config[job2idx:]:
            f2.write(line + '\n')

    with open(subfile_1, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch_1, command_1, copy_the_refine_xyz_1, clean_scratch))

    with open(subfile_2, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch_2, command_2, copy_the_refine_xyz_2, clean_scratch))

    return subfile_1, subfile_2


def update_irc_opt_status(qm_collection:object, target:object, job_id_1:str, job_id_2:str):
    update_field = {'irc_forward_opt_status': "job_launched", 'irc_forward_opt_jobid': job_id_1,
                    'irc_backward_opt_status': "job_launched", 'irc_backward_opt_jobid': job_id_2}
    qm_collection.update_one(target, {"$unset": {'irc_opt_status': ""}, "$set": update_field}, True)

"""
Submmit QMMM calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def launch_qmmm_opt_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_opt')
    count = 0
    for target in targets[:num]:
        qmmm_reactant_dir = path.join(target['path'], 'QMMM_REACTANT/')
        IRC_dir_path = path.join(target['path'], 'IRC/')
        if target['irc_equal'] == 'forward equal to reactant':
            reactant = path.join(IRC_dir_path, 'irc_forward.xyz')
            product = path.join(IRC_dir_path, 'irc_backward.xyz')
        else:
            reactant = path.join(IRC_dir_path, 'irc_backward.xyz')
            product = path.join(IRC_dir_path, 'irc_forward.xyz')
        if os.path.exists(qmmm_reactant_dir):
            os.chdir(qmmm_reactant_dir)
        else:
            os.mkdir(qmmm_reactant_dir)
            os.chdir(qmmm_reactant_dir)
        subfile_1 = create_qmmm_opt(qmmm_reactant_dir, reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd_1 = 'qsub {}'.format(subfile_1)
        process = subprocess.Popen([cmd_1],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_1 = stdout.decode().replace("\n", "")

        qmmm_product_dir = path.join(target['path'], 'QMMM_PRODUCT/')
        if os.path.exists(qmmm_product_dir):
            os.chdir(qmmm_product_dir)
        else:
            os.mkdir(qmmm_product_dir)
            os.chdir(qmmm_product_dir)
        subfile_2 = create_qmmm_opt(qmmm_product_dir, product, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd_2 = 'qsub {}'.format(subfile_2)
        process = subprocess.Popen([cmd_2],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_2 = stdout.decode().replace("\n", "")
        # update status job_launched
        update_qmmm_opt_status(qm_collection, target, job_id_1, job_id_2)
        count += 1
    print(highlight_text('QMMM OPT'))
    print('\nQMMM opt launced {} jobs (forward + backward)\n'.format(count * 2))

def create_qmmm_opt(qmmm_dir:str, target_geometry:str, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(qmmm_dir)))), 'config')
    qmmm_opt_config = path.join(base_dir_path, 'qmmm_opt.lot')
    qmmm_opt_input = path.join(qmmm_dir, 'qmmm_opt.in')
    qmmm_opt_output = path.join(qmmm_dir, 'qmmm_opt.out')
    subfile = path.join(qmmm_dir, 'qmmm_opt.job')

    with open(target_geometry, 'r') as f1:
        lines = f1.read().splitlines()
    with open(qmmm_opt_config) as f:
        config = [line.strip() for line in f]
    qm_atoms = []
    qm_xyzs = []
    for i, text in enumerate(config):
        if text.upper().startswith('$QM_ATOMS'):
            break
    for qm_atom in config[i+1:]:
        if '$end' in qm_atom:
            break
        else:
            try:  
                qm_atoms.append(int(qm_atom))
            except:
                continue
    for j, text in enumerate(config):
        if text.upper().startswith('$MOLECULE'):
            break
    nqm_atoms = len(qm_atoms)
    for qm_xyz in config[j+1: j+nqm_atoms+2]:
        qm_xyz = qm_xyz.split()
        qm_xyzs.append(qm_xyz)
    with open(qmmm_opt_input, 'w') as f:
        for k, text in enumerate(config):
            if '$MOLECULE' not in text.upper():
                f.write(text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for l, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+l] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for last_text in config[k + 2 + nqm_atoms:]:
            f.write(last_text + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(qmmm_dir)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, qmmm_opt_input, qmmm_opt_output)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))
    return subfile

def update_qmmm_opt_status(qm_collection:object, target:object, job_id_1:str, job_id_2:str):
    update_field = {'qmmm_opt_reactant_status': "job_launched", 'qmmm_opt_reactant_jobid': job_id_1,
                    'qmmm_opt_product_status': "job_launched", 'qmmm_opt_product_jobid': job_id_2}
    qm_collection.update_one(target, {"$unset": {'qmmm_opt_status': ""}, "$set": update_field}, True)

"""
QMMM FREQ OPT
"""

def launch_qmmm_freq_opt_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_freq_opt')
    count = 0
    for target in targets[:num]:
        qmmm_reactant_dir = path.join(target['path'], 'QMMM_REACTANT/')
        reactant = path.join(qmmm_reactant_dir, 'qmmm_opt.xyz')
        if os.path.exists(qmmm_reactant_dir):
            os.chdir(qmmm_reactant_dir)
        else:
            os.mkdir(qmmm_reactant_dir)
            os.chdir(qmmm_reactant_dir)
        subfile_1 = create_qmmm_freq_opt(qmmm_reactant_dir, reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd_1 = 'qsub {}'.format(subfile_1)
        process = subprocess.Popen([cmd_1],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_1 = stdout.decode().replace("\n", "")

        qmmm_product_dir = path.join(target['path'], 'QMMM_PRODUCT/')
        product = path.join(qmmm_product_dir, 'qmmm_opt.xyz')
        if os.path.exists(qmmm_product_dir):
            os.chdir(qmmm_product_dir)
        else:
            os.mkdir(qmmm_product_dir)
            os.chdir(qmmm_product_dir)
        subfile_2 = create_qmmm_freq_opt(qmmm_product_dir, product, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd_2 = 'qsub {}'.format(subfile_2)
        process = subprocess.Popen([cmd_2],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_2 = stdout.decode().replace("\n", "")
        # update status job_launched
        update_qmmm_freq_opt_status(qm_collection, target, job_id_1, job_id_2)
        count += 1
    print(highlight_text('QMMM FREQ OPT'))
    print('\nQMMM freq opt launced {} jobs (forward + backward)\n'.format(count * 2))

def launch_qmmm_freq_opt_restart_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_freq_opt_reactant')
    count = 0
    for target in targets[:num]:
        if target['qmmm_freq_opt_reactant_status'] == 'restart':
            qmmm_reactant_dir = path.join(target['path'], 'QMMM_REACTANT/')
            reactant = path.join(qmmm_reactant_dir, 'qmmm_freq_opt.xyz')
            subfile_1 = create_qmmm_freq_opt(qmmm_reactant_dir, reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
            cmd_1 = 'qsub {}'.format(subfile_1)
            process = subprocess.Popen([cmd_1],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, shell=True)
            stdout, stderr = process.communicate()
            # get job id from stdout, e.g., "106849.h81"
            job_id_1 = stdout.decode().replace("\n", "")
            update_status(qm_collection, target, job_name='qmmm_freq_opt_reactant', job_id=job_id_1)
            count += 1

    targets = select_targets(qm_collection, job_name='qmmm_freq_opt_product')

    for target in targets[:num]:
        if target['qmmm_freq_opt_product_status'] == 'restart':
            qmmm_product_dir = path.join(target['path'], 'QMMM_PRODUCT/')
            product = path.join(qmmm_product_dir, 'qmmm_freq_opt.xyz')
            subfile_2 = create_qmmm_freq_opt(qmmm_product_dir, product, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
            cmd_2 = 'qsub {}'.format(subfile_2)
            process = subprocess.Popen([cmd_1],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, shell=True)
            stdout, stderr = process.communicate()
            # get job id from stdout, e.g., "106849.h81"
            job_id_2 = stdout.decode().replace("\n", "")
            update_status(qm_collection, target, job_name='qmmm_freq_opt_product', job_id=job_id_2)
            count += 1

    print(highlight_text('QMMM FREQ OPT RESTART'))
    print('\nQMMM freq opt restart launced {} jobs\n'.format(count))

def create_qmmm_freq_opt(qmmm_dir:str, target_geometry:str, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(qmmm_dir)))), 'config')
    qmmm_freq_opt_config = path.join(base_dir_path, 'qmmm_freq_opt.lot')
    qmmm_freq_opt_input = path.join(qmmm_dir, 'qmmm_freq_opt.in')
    qmmm_freq_opt_output = path.join(qmmm_dir, 'qmmm_freq_opt.out')
    subfile = path.join(qmmm_dir, 'qmmm_freq_ts.job')

    with open(target_geometry, 'r') as f1:
        lines = f1.read().splitlines()
    with open(qmmm_freq_opt_config) as f:
        config = [line.strip() for line in f]
    qm_atoms = []
    qm_xyzs = []
    for i, text in enumerate(config):
        if text.upper().startswith('$QM_ATOMS'):
            break
    for qm_atom in config[i+1:]:
        if '$end' in qm_atom:
            break
        else:
            try:  
                qm_atoms.append(int(qm_atom))
            except:
                continue
    for j, text in enumerate(config):
        if text.upper().startswith('$MOLECULE'):
            break
    nqm_atoms = len(qm_atoms)
    for qm_xyz in config[j+1: j+nqm_atoms+2]:
        qm_xyz = qm_xyz.split()
        qm_xyzs.append(qm_xyz)
    with open(qmmm_freq_opt_input, 'w') as f:
        for k, text in enumerate(config):
            if '$MOLECULE' not in text.upper():
                f.write(text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for l, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+l] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for m, last_text in enumerate(config[k + 2 + nqm_atoms:]):
            if '$MOLECULE' not in last_text.upper():
                f.write(last_text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for n, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+n] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for last_text in config[k + m + 2 * 2 + nqm_atoms * 2:]:
            f.write(last_text + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(qmmm_dir)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, qmmm_freq_opt_input, qmmm_freq_opt_output)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))
    
    return subfile

def update_qmmm_freq_opt_status(qm_collection:object, target:object, job_id_1:str, job_id_2:str):
    update_field = {'qmmm_freq_opt_reactant_status': "job_launched", 'qmmm_freq_opt_reactant_jobid': job_id_1,
                    'qmmm_freq_opt_product_status': "job_launched", 'qmmm_freq_opt_product_jobid': job_id_2}
    qm_collection.update_one(target, {"$unset": {'qmmm_opt_status': ""}, "$set": update_field}, True)

"""
QMMM FREQ TS
"""

def launch_qmmm_freq_ts_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_freq_ts')
    count = 0
    for target in targets[:num]:
        qmmm_ts_dir = path.join(target['path'], 'QMMM_TS/')
        IRC_dir_path = path.join(target['path'], 'IRC/')
        if target['qmmm_freq_ts_status'] == 'restart':
            ts = path.join(qmmm_ts_dir, 'qmmm_ts.xyz')
        else:
            ts = path.join(IRC_dir_path, 'ts_geo.xyz')

        if os.path.exists(qmmm_ts_dir):
            os.chdir(qmmm_ts_dir)
        else:
            os.mkdir(qmmm_ts_dir)
            os.chdir(qmmm_ts_dir)

        subfile = create_qmmm_freq_ts(qmmm_ts_dir, ts, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        update_status(qm_collection, target, job_name = 'qmmm_freq_ts', job_id = job_id)
        count += 1
    print(highlight_text('QMMM TS'))
    print('\nQMMM ts launced {} jobs\n'.format(count))

def create_qmmm_freq_ts(qmmm_dir:str, target_geometry:str, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(qmmm_dir)))), 'config')
    qmmm_opt_config = path.join(base_dir_path, 'qmmm_freq_ts.lot')
    qmmm_freq_ts_input = path.join(qmmm_dir, 'qmmm_freq_ts.in')
    qmmm_freq_ts_output = path.join(qmmm_dir, 'qmmm_freq_ts.out')
    subfile = path.join(qmmm_dir, 'qmmm_freq_ts.job')

    with open(target_geometry, 'r') as f1:
        lines = f1.read().splitlines()
    with open(qmmm_opt_config) as f:
        config = [line.strip() for line in f]
    qm_atoms = []
    qm_xyzs = []
    for i, text in enumerate(config):
        if text.upper().startswith('$QM_ATOMS'):
            break
    for qm_atom in config[i+1:]:
        if '$end' in qm_atom:
            break
        else:
            try:  
                qm_atoms.append(int(qm_atom))
            except:
                continue
    for j, text in enumerate(config):
        if text.upper().startswith('$MOLECULE'):
            break
    nqm_atoms = len(qm_atoms)
    for qm_xyz in config[j+1: j+nqm_atoms+2]:
        qm_xyz = qm_xyz.split()
        qm_xyzs.append(qm_xyz)
    with open(qmmm_freq_ts_input, 'w') as f:
        for k, text in enumerate(config):
            if '$MOLECULE' not in text.upper():
                f.write(text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for l, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+l] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for m, last_text in enumerate(config[k + 2 + nqm_atoms:]):
            if '$MOLECULE' not in last_text.upper():
                f.write(last_text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for n, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+n] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for last_text in config[k + m + 2 * 2 + nqm_atoms * 2:]:
            f.write(last_text + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(qmmm_dir)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, qmmm_freq_ts_input, qmmm_freq_ts_output)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))

    return subfile

"""
QMMM FREQ
After QMMM opt then run a freq to check if exist imaginary frequency
"""

def launch_qmmm_freq_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_freq')
    count = 0
    for target in targets[:num]:
        qmmm_reactant_dir = path.join(target['path'], 'QMMM_REACTANT/')
        reactant = path.join(qmmm_reactant_dir, 'qmmm_freq_opt.xyz')

        if os.path.exists(qmmm_reactant_dir):
            os.chdir(qmmm_reactant_dir)
        else:
            os.mkdir(qmmm_reactant_dir)
            os.chdir(qmmm_reactant_dir)

        subfile_1 = create_qmmm_freq(qmmm_reactant_dir, reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd_1 = 'qsub {}'.format(subfile_1)
        process = subprocess.Popen([cmd_1],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_1 = stdout.decode().replace("\n", "")

        qmmm_product_dir = path.join(target['path'], 'QMMM_PRODUCT/')
        product = path.join(qmmm_product_dir, 'qmmm_freq_opt.xyz')
        if os.path.exists(qmmm_product_dir):
            os.chdir(qmmm_product_dir)
        else:
            os.mkdir(qmmm_product_dir)
            os.chdir(qmmm_product_dir)
        subfile_2 = create_qmmm_freq(qmmm_product_dir, product, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd_2 = 'qsub {}'.format(subfile_2)
        process = subprocess.Popen([cmd_2],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_2 = stdout.decode().replace("\n", "")
        # update status job_launched
        update_qmmm_freq_status(qm_collection, target, job_id_1, job_id_2)
        count += 1
    print(highlight_text('QMMM FREQ'))
    print('\nQMMM freq launced {} jobs (forward + backward)\n'.format(count * 2))

def create_qmmm_freq(qmmm_dir:str, target_geometry:str, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(qmmm_dir)))), 'config')
    qmmm_freq_config = path.join(base_dir_path, 'qmmm_freq.lot')
    qmmm_freq_input = path.join(qmmm_dir, 'qmmm_freq.in')
    qmmm_freq_output = path.join(qmmm_dir, 'qmmm_freq.out')
    subfile = path.join(qmmm_dir, 'qmmm_freq.job')

    with open(target_geometry, 'r') as f1:
        lines = f1.read().splitlines()
    with open(qmmm_freq_config) as f:
        config = [line.strip() for line in f]
    qm_atoms = []
    qm_xyzs = []
    for i, text in enumerate(config):
        if text.upper().startswith('$QM_ATOMS'):
            break
    for qm_atom in config[i+1:]:
        if '$end' in qm_atom:
            break
        else:
            try:  
                qm_atoms.append(int(qm_atom))
            except:
                continue
    for j, text in enumerate(config):
        if text.upper().startswith('$MOLECULE'):
            break
    nqm_atoms = len(qm_atoms)
    for qm_xyz in config[j+1: j+nqm_atoms+2]:
        qm_xyz = qm_xyz.split()
        qm_xyzs.append(qm_xyz)
    with open(qmmm_freq_input, 'w') as f:
        for k, text in enumerate(config):
            if '$MOLECULE' not in text.upper():
                f.write(text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for l, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+l] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for last_text in config[k + 2 + nqm_atoms:]:
            f.write(last_text + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(qmmm_dir)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, qmmm_freq_input, qmmm_freq_output)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))

    return subfile

def update_qmmm_freq_status(qm_collection:object, target:object, job_id_1:str, job_id_2:str):
    update_field = {'qmmm_freq_reactant_status': "job_launched", 'qmmm_freq_reactant_jobid': job_id_1,
                    'qmmm_freq_product_status': "job_launched", 'qmmm_freq_product_jobid': job_id_2}
    qm_collection.update_one(target, {"$unset": {'qmmm_freq_status': ""}, "$set": update_field}, True)

"""
QMMM TS FREQ
After QMMM ts converge then run a freq
"""

def launch_qmmm_ts_freq_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_targets(qm_collection, job_name='qmmm_ts_freq')
    count = 0
    for target in targets[:num]:
        qmmm_ts_dir = path.join(target['path'], 'QMMM_TS/')
        ts = path.join(qmmm_ts_dir, 'qmmm_ts.xyz')

        if os.path.exists(qmmm_ts_dir):
            os.chdir(qmmm_ts_dir)
        else:
            os.mkdir(qmmm_ts_dir)
            os.chdir(qmmm_ts_dir)
        
        subfile = create_qmmm_freq(qmmm_ts_dir, ts, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        update_status(qm_collection, target, job_name = 'qmmm_ts_freq', job_id = job_id)
        count += 1
    print(highlight_text('QMMM TS FREQ'))
    print('\nQMMM ts freq launced {} jobs\n'.format(count))

"""
QMMM REFINE
"""

def select_qmmm_refine_targets(qm_collection:object) -> list:
    query = {'$and':
            [{"qmmm_refine_status":
                {"$in":
                    ["job_unrun"]}
                },
                {"qmmm_ts_refine":
                {"$in":
                    ["job_unrun"]}
                }]}
    targets = list(qm_collection.find(query))
    return targets

def launch_qmmm_refine_jobs(qm_collection:object, num:int=10, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16):
    targets = select_qmmm_refine_targets(qm_collection)
    count = 0
    for target in targets[:num]:
        qmmm_reactant_dir = path.join(target['path'], 'QMMM_REACTANT/')
        reactant = path.join(qmmm_reactant_dir, 'qmmm_final.xyz')
        if os.path.exists(qmmm_reactant_dir):
            os.chdir(qmmm_reactant_dir)
        else:
            os.mkdir(qmmm_reactant_dir)
            os.chdir(qmmm_reactant_dir)
        subfile_1 = create_qmmm_refine(qmmm_reactant_dir, reactant, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        
        cmd_1 = 'qsub {}'.format(subfile_1)
        process = subprocess.Popen([cmd_1],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_1 = stdout.decode().replace("\n", "")

        qmmm_product_dir = path.join(target['path'], 'QMMM_PRODUCT/')
        product = path.join(qmmm_product_dir, 'qmmm_final.xyz')
        if os.path.exists(qmmm_product_dir):
            os.chdir(qmmm_product_dir)
        else:
            os.mkdir(qmmm_product_dir)
            os.chdir(qmmm_product_dir)
        subfile_2 = create_qmmm_refine(qmmm_product_dir, product, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd_2 = 'qsub {}'.format(subfile_2)
        process = subprocess.Popen([cmd_2],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_2 = stdout.decode().replace("\n", "")

        qmmm_ts_dir = path.join(target['path'], 'QMMM_TS/')
        ts = path.join(qmmm_ts_dir, 'qmmm_final.xyz')
        if os.path.exists(qmmm_ts_dir):
            os.chdir(qmmm_ts_dir)
        else:
            os.mkdir(qmmm_ts_dir)
            os.chdir(qmmm_ts_dir)
        subfile_3 = create_qmmm_refine(qmmm_ts_dir, ts, ncpus=ncpus, mpiprocs=mpiprocs, ompthreads=ompthreads)
        cmd_3 = 'qsub {}'.format(subfile_3)
        process = subprocess.Popen([cmd_3],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id_3 = stdout.decode().replace("\n", "")

        # update status job_launched
        update_qmmm_refine_status(qm_collection, target, job_id_1, job_id_2, job_id_3)
        count += 1
    print(highlight_text('QMMM REFINE'))
    print('\nQMMM refine launced {} jobs (reactant + product + ts)\n'.format(count * 3))

def create_qmmm_refine(qmmm_dir:str, target_geometry:str, ncpus:int=16, mpiprocs:int=1, ompthreads:int=16) -> str:
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(qmmm_dir)))), 'config')
    qmmm_sp_config = path.join(base_dir_path, 'qmmm_sp.lot')
    qmmm_sp_input = path.join(qmmm_dir, 'qmmm_sp.in')
    qmmm_sp_output = path.join(qmmm_dir, 'qmmm_sp.out')
    subfile = path.join(qmmm_dir, 'qmmm_sp.job')

    with open(target_geometry, 'r') as f1:
        lines = f1.read().splitlines()
    with open(qmmm_sp_config) as f:
        config = [line.strip() for line in f]
    qm_atoms = []
    qm_xyzs = []
    for i, text in enumerate(config):
        if text.upper().startswith('$QM_ATOMS'):
            break
    for qm_atom in config[i+1:]:
        if '$end' in qm_atom:
            break
        else:
            try:  
                qm_atoms.append(int(qm_atom))
            except:
                continue
    for j, text in enumerate(config):
        if text.upper().startswith('$MOLECULE'):
            break
    nqm_atoms = len(qm_atoms)
    for qm_xyz in config[j+1: j+nqm_atoms+2]:
        qm_xyz = qm_xyz.split()
        qm_xyzs.append(qm_xyz)
    with open(qmmm_sp_input, 'w') as f:
        for k, text in enumerate(config):
            if '$MOLECULE' not in text.upper():
                f.write(text + '\n')
            else:
                break
        f.write('$MOLECULE' + '\n')
        for l, line in enumerate(qm_xyzs):
            if len(line) > 2:
                connectivity = '\t'.join(line[-5:])
                geometry = lines[1+l] + ' \t ' + connectivity
                f.write(geometry + '\n')
            else:
                line = ' '.join(line)
                f.write(line + '\n')
        for last_text in config[k + 2 + nqm_atoms:]:
            f.write(last_text + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nstart=$(date +\'%s\')\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(qmmm_dir)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, qmmm_sp_input, qmmm_sp_output)
    clean_scratch = """rm -r $QCSCRATCH\necho "It took $(($(date +'%s') - $start)) seconds" """

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))
    return subfile

def update_qmmm_refine_status(qm_collection:object, target:object, job_id_1:str, job_id_2:str, job_id_3:str):
    update_field = {'qmmm_refine_reactant_status': "job_launched", 'qmmm_refine_reactant_jobid': job_id_1,
                    'qmmm_refine_product_status': "job_launched", 'qmmm_refine_product_jobid': job_id_2,
                    'qmmm_refine_ts_status': "job_launched", 'qmmm_refine_ts_jobid': job_id_3}
    qm_collection.update_one(target, {"$unset": {'qmmm_refine_status': "", 'qmmm_ts_refine':""}, "$set": update_field}, True)


def launch_jobs(num=30, level_of_theory='ORCA', ncpus=4, mpiprocs=1, ompthreads=4):
    # Hcap = {} mean the last {} H atoms
    qm_collection = db['qm_calculate_center']
    launch_ssm_jobs(qm_collection, num=50, level_of_theory='ORCA',ncpus=1, mpiprocs=1, ompthreads=1)
    launch_ts_refine_jobs(qm_collection, num=50, ncpus=1, mpiprocs=1, ompthreads=1)
    launch_ts_jobs(qm_collection, num=10, level_of_theory=level_of_theory, ncpus=4, mpiprocs=4, ompthreads=1, Hcap=12)
    launch_irc_jobs(qm_collection, num=10, ncpus=8, mpiprocs=8, ompthreads=1)
    launch_irc_opt_jobs(qm_collection, num=10, level_of_theory=level_of_theory,ncpus=8, mpiprocs=8, ompthreads=1, Hcap=12)

    # launch_qmmm_opt_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_freq_opt_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_freq_opt_restart_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_freq_ts_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_freq_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_ts_freq_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)
    # launch_qmmm_refine_jobs(qm_collection, num=num, ncpus=16, mpiprocs=1, ompthreads=16)

print_header()
launch_jobs(num=3, level_of_theory='ORCA', ncpus=4, mpiprocs=4, ompthreads=1)
