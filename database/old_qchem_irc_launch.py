"""
####
IRC in qchem is not that robust...
So, we switch to pysisyphus
The following may be droped or just remain for someday qchem's irc work well.
Otherwise, make an interface between pysisyphus and qchem.
Now pysisyphus only support gaussian, orca, xtb, etc......
####


Submmit IRC(Intrinsic Reaction Coordinate, in Qchem called 'rpath') calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"

    
def select_irc_target():
    qm_collection = db['qm_calculate_center']
    reg_query = {"irc_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_irc_jobs():
    targets = select_irc_target()
    
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.mkdir(IRC_dir_path)
        os.chdir(IRC_dir_path)
            
        TS_dir_path = path.join(target, 'TS/')
        subfile_1,  subfile_2= create_irc_sub_file(TS_dir_path, IRC_dir_path)
        cmd_1 = 'qsub {}'.format(subfile_1)
        process_1 = subprocess.Popen([cmd_1],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process_1.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_status(target, job_id, direction = 'forward')
        
        cmd_2 = 'qsub {}'.format(subfile_2)
        process_2 = subprocess.Popen([cmd_2],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process_2.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_status(target, job_id, direction = 'reverse')
        
def create_irc_sub_file(TS_dir_path, IRC_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    ts_geo_path = path.join(TS_dir_path, 'ts_geo.xyz')
    
    irc_forward_input_file = path.join(IRC_dir_path, 'irc_forward.in')
    irc_reverse_input_file = path.join(IRC_dir_path, 'irc_reverse.in')
    irc_forward_output_file = path.join(IRC_dir_path, 'irc_forward.out')
    irc_reverse_output_file = path.join(IRC_dir_path, 'irc_reverse.out')

    subfile_1 = path.join(IRC_dir_path, 'cal_irc_forward.job')
    subfile_2 = path.join(IRC_dir_path, 'cal_irc_reverse.job')

    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(IRC_dir_path)))), 'config')
    irc_forward_lot = path.join(base_dir_path, 'freq_irc_forward.lot')
    irc_reverse_lot = path.join(base_dir_path, 'freq_irc_reverse.lot')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(IRC_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command_1 = 'qchem -nt {} {} {}'.format(ncpus, irc_forward_input_file, irc_forward_output_file)
    command_2 = 'qchem -nt {} {} {}'.format(ncpus, irc_reverse_input_file, irc_reverse_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(irc_forward_lot) as f:
        forward_config = [line.strip() for line in f]
    with open(irc_reverse_lot) as f:
        reverse_config = [line.strip() for line in f]

    with open(ts_geo_path, 'r') as f1:
        lines = f1.read().splitlines()

    with open(irc_forward_input_file, 'w') as f2:
        for i, text in enumerate(forward_config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                forward_config[(i+1):(i+1)] = cblock
                break
        for line in forward_config:
            f2.write(line + '\n')
    with open(irc_reverse_input_file, 'w') as f2:
        for i, text in enumerate(reverse_config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                reverse_config[(i+1):(i+1)] = cblock
                break
        for line in reverse_config:
            f2.write(line + '\n')

    with open(subfile_1, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command_1, clean_scratch))
    with open(subfile_2, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command_2, clean_scratch))
            
    return subfile_1, subfile_2
    

def update_irc_status(target, job_id, direction):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    irc_status = 'irc_{}_status'.format(str(direction))
    irc_jobid = 'irc_{}_jobid'.format(str(direction))
    update_field = {irc_status:"job_launched", irc_jobid:job_id}
    qm_collection.update_one(reg_query, {"$unset": {'irc_status':""}, "$set": update_field}, True)


Submmit opt job which is from irc
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"


def select_irc_opt_target(direction = 'forward'):
    
    qm_collection = db['qm_calculate_center']
    irc_opt_status = 'opt_{}_status'.format(direction)
    reg_query = {irc_opt_status:"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_irc_opt_jobs():
    
    targets = select_irc_opt_target(direction = 'forward')
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.chdir(IRC_dir_path)
        subfile = create_irc_opt_sub_file(IRC_dir_path, direction = 'forward')
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(target, job_id, direction = 'forward')
        
    targets = select_irc_opt_target(direction = 'reverse')
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.chdir(IRC_dir_path)
        subfile = create_irc_opt_sub_file(IRC_dir_path, direction = 'reverse')
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(target, job_id, direction = 'reverse')


def update_irc_opt_status(target, job_id, direction):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    irc_opt_status = 'opt_{}_status'.format(direction)
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str(direction))
    update_field = {irc_opt_status:"opt_job_launched", irc_opt_jobid:job_id}
    qm_collection.update_one(reg_query, {"$set": update_field}, True)
    

def create_irc_opt_sub_file(irc_path, direction = 'forward', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    job_name = 'irc_{}_opt.job'.format(direction)
    subfile = path.join(irc_path, job_name)
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(irc_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    inputname = '{}_opt.in'.format(direction)
    outputname = '{}_opt.out'.format(direction)
    nes4 = 'qchem -nt {} {} {}'.format(ncpus, inputname, outputname)
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
    return subfile
"""
