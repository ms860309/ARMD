"""
IRC check.
"""


def select_irc_target(direction='forward'):
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    irc_status = 'irc_{}_status'.format(direction)
    query = {irc_status:
             {"$in":
              ["job_launched", "job_running", "job_queueing"]
              }
             }
    targets = list(qm_collection.find(query))

    return targets


def check_irc_status(job_id):
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
        return "job_queueing"
    else:
        return "job_launched"


def check_irc_content(target_path, direction='forward'):
    reactant_path = path.join(target_path, 'reactant.xyz')
    irc_path = path.join(target_path, 'IRC/')
    base_dir_path = path.join(path.dirname(path.dirname(
        path.dirname(path.dirname(irc_path)))), 'config')
    opt_lot = path.join(base_dir_path, 'opt_freq.lot')
    opt_name = '{}_opt.in'.format(direction)
    opt_in = path.join(irc_path, opt_name)
    irc_output = path.join(irc_path, 'irc_{}.out'.format(direction))

    with open(opt_lot) as f:
        config = [line.strip() for line in f]
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    with open(irc_output, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 12280, 0), 0)  # Read last 12 kB of file
        lines = f.readlines()

    if lines[-2] == ' IRC backup failure\n' or lines[-2] == ' IRC failed final bisector step\n':

        with open(irc_output, 'r') as f:
            full_lines = f.readlines()

        # Sometimes irc success in forward(reverse) but fail in reverse(forward).
        # We wan't to catch the final structure to optimize.
        # But if "convergence criterion reached" in first direction fail in second that will cause reactant equal to product.

        if '  IRC -- convergence criterion reached.\n' in full_lines:
            for idx, i in enumerate(full_lines):
                if i.startswith('  IRC -- convergence criterion reached.\n'):
                    break
            full_lines = full_lines[:idx]
            for idx2, j in enumerate(reversed(full_lines)):
                if j.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
                    break
            geo = []
            for i in full_lines[-idx2 + 2: -idx2 + 2 + atom_number]:
                atom = i.split()[1:]
                geo.append('  '.join(atom))
            with open(opt_in, 'w') as f:
                f.write('$molecule\n{} {}\n'.format(0, 1))
                f.write('\n'.join(geo))
                f.write('\n$end\n\n')
                for line in config:
                    f.write(line + '\n')
        else:
            for idx, i in enumerate(lines):
                if i.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
                    break
            geo = []
            for i in lines[idx + 3: idx + 3 + atom_number]:
                atom = i.split()[1:]
                geo.append('  '.join(atom))

            with open(opt_in, 'w') as f:
                f.write('$molecule\n{} {}\n'.format(0, 1))
                f.write('\n'.join(geo))
                f.write('\n$end\n\n')
                for line in config:
                    f.write(line + '\n')

        return 'need opt'
    elif lines[-5] == '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
        return 'job_success'
    elif lines[-2] == ' Bad initial gradient\n':
        return 'Bad initial gradient'
    elif lines[-2] == ' IRC --- Failed line search\n':
        return 'Failed line search'
    elif lines[-6] == ' Error in gen_scfman\n':
        return 'Error in gen_scfman'
    else:
        return 'unknown fail information'


def generate_irc_product_xyz(target, direction='forward'):
    irc_path = path.join(target['path'], 'IRC')
    reactant_path = path.join(target['path'], 'reactant.xyz')
    output_name = 'irc_{}.out'.format(direction)
    output = path.join(irc_path, output_name)
    name = path.join(irc_path, '{}.xyz'.format(direction))

    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])

    with open(output, 'r') as f:
        full_lines = f.readlines()
    count = 1
    for idx, i in enumerate(full_lines):
        if i.startswith('  IRC -- convergence criterion reached.\n'):
            count += 1
            if count == 2:
                break
    full_lines = full_lines[:idx]
    for idx2, j in enumerate(reversed(full_lines)):
        if j.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
            break
    geo = []
    for i in full_lines[-idx2 + 2: -idx2 + 2 + atom_number]:
        atom = i.split()[1:]
        geo.append('  '.join(atom))
    with open(name, 'w') as f:
        f.write(str(atom_number))
        f.write('\n\n')
        f.write('\n'.join(geo))


def check_irc_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_irc_target(direction='forward')

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['irc_forward_jobid']
        # 2. check the job pbs status
        new_status = check_irc_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_content(target['path'], direction='forward')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('forward')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_irc_product_xyz(target, direction='forward')
                update_field = {
                    irc_status: new_status, 'irc_equal': 'waiting for check'
                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('forward')
                update_field = {
                    irc_status: new_status, opt_status: 'job_unrun'
                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                    irc_status: new_status
                }
                qm_collection.update_one(target, {"$set": update_field}, True)

    # 1. select jobs to check
    targets = select_irc_target(direction='reverse')

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['irc_reverse_jobid']
        # 2. check the job pbs status
        new_status = check_irc_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_content(target['path'], direction='reverse')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('reverse')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_irc_product_xyz(target, direction='reverse')
                update_field = {
                    irc_status: new_status, 'irc_equal': 'waiting for check'
                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('reverse')
                update_field = {
                    irc_status: new_status, opt_status: 'job_unrun'
                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                    irc_status: new_status
                }
                qm_collection.update_one(target, {"$set": update_field}, True)


"""
IRC  success check.
This is to check whether irc forward and reverse direction equal to expectation.
"""


def select_irc_equal_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    query = {'$and':
             [
                 {"irc_forward_status":
                  {"$in":
                   ['job_success', 'opt_success']}},
                 {'irc_reverse_status':
                  {'$in':
                   ['job_success', 'opt_success']}},
                 {'irc_equal':
                  {'$in':
                   ['waiting for check']}}
             ]
             }
    targets = list(qm_collection.find(query))

    return targets


def check_irc_equal():

    targets = select_irc_equal_target()
    qm_collection = db['qm_calculate_center']

    for target in targets:
        if target['irc_forward_status'] in ['job_success', 'opt_success'] and target['irc_reverse_status'] in ['job_success', 'opt_success']:
            new_status = check_irc_equal_status(target)
            orig_status = target['irc_equal']
            if orig_status != new_status:
                if new_status == 'forward equal to reactant and reverse equal to product' or new_status == 'reverse equal to reactant and forward equal to product' or new_status == 'reverse equal to reactant but forward does not equal to product' or new_status == 'forward equal to reactant but reverse does not equal to product':
                    update_field = {
                        'irc_equal': new_status, 'energy_status': 'job_unrun'
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)
                else:
                    update_field = {
                        'irc_equal': new_status
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)
        elif target['opt_reverse_status'] == 'job_fail' or target['opt_forward_status'] == 'job_fail':
            update_field = {
                'irc_equal': 'opt fail'
            }
            qm_collection.update_one(target, {"$set": update_field}, True)


def check_irc_equal_status(target):

    irc_path = path.join(target['path'], 'IRC/')
    reactant_path = path.join(target['path'], 'reactant.xyz')
    product_path = path.join(target['path'], 'ssm_product.xyz')
    forward_output = path.join(irc_path, 'forward.xyz')
    reverse_output = path.join(irc_path, 'reverse.xyz')

    pyMol_1 = xyz_to_pyMol(reactant_path)
    pyMol_2 = xyz_to_pyMol(product_path)
    pyMol_3 = xyz_to_pyMol(forward_output)
    pyMol_4 = xyz_to_pyMol(reverse_output)

    if pyMol_3.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward equal to reverse'
    elif (pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip()) and (pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip()):
        return 'forward and reverse equal to reactant'
    elif (pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip()) and (pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip()):
        return 'forward and reverse equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward equal to reactant and reverse equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'reverse equal to reactant and forward equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() != pyMol_3.write('inchiKey').strip():
        return 'reverse equal to reactant but forward does not equal to product'
    elif pyMol_1.write('inchiKey').strip() != pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'reverse does not equal to reactant but forward equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() != pyMol_4.write('inchiKey').strip():
        return 'forward equal to reactant but reverse does not equal to product'
    elif pyMol_1.write('inchiKey').strip() != pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward does not equal to reactant but reverse equal to product'
    else:
        return 'unknown (Maybe both of them are not equal to reactant&product)'


"""
IRC opt check.
"""


def select_irc_opt_target(direction='forward'):
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    irc_opt_status = 'opt_{}_status'.format(direction)
    reg_query = {irc_opt_status:
                 {"$in":
                  ["opt_job_launched", "opt_job_running", "opt_job_queueing"]
                  }
                 }
    targets = list(qm_collection.find(reg_query))

    return targets


def check_irc_opt_job_status(job_id):
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
        return "opt_job_running"
    elif stdout[idx+2] == 'Q':
        return 'opt_job_queueing'
    else:
        return "opt_job_launched"


def check_irc_opt_content(dir_path, direction='forward'):
    reactant_path = os.path.join(dir_path, 'reactant.xyz')
    irc_path = path.join(dir_path, "IRC")
    xyzname = '{}.xyz'.format(direction)
    output = path.join(irc_path, xyzname)
    output_path = path.join(irc_path, '{}_opt.out'.format(direction))

    try:
        q = QChem(outputfile=output_path)
        q.create_geo_file(output)
        return 'job_success'
    except:
        return 'job_fail'


def check_irc_opt_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_irc_opt_target(direction='forward')
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str('forward'))
    qm_collection = db['qm_calculate_center']

    # 2. check the job pbs_status
    for target in targets:
        job_id = target[irc_opt_jobid]
        new_status = check_irc_opt_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_opt_content(
                target['path'], direction='forward')

            # 4. check with original status which
            # should be job_launched or job_running
            # if any difference update status
            irc_status = 'irc_{}_status'.format(str('forward'))
            irc_opt_status = 'opt_{}_status'.format(str('forward'))
            orig_status = target[irc_opt_status]
            if orig_status != new_status:
                if new_status == 'job_success':
                    update_field = {
                        irc_status: 'opt_success', irc_opt_status: new_status, 'irc_equal': 'waiting for check'
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)
                else:
                    update_field = {
                        irc_status: 'opt_fail', irc_opt_status: new_status
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)

    # 1. select jobs to check
    targets = select_irc_opt_target(direction='reverse')
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str('reverse'))

    # 2. check the job pbs_status
    for target in targets:
        job_id = target[irc_opt_jobid]
        new_status = check_irc_opt_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_opt_content(
                target['path'], direction='reverse')

            # 4. check with original status which
            # should be job_launched or job_running
            # if any difference update status
            irc_status = 'irc_{}_status'.format(str('reverse'))
            irc_opt_status = 'opt_{}_status'.format(str('reverse'))
            orig_status = target[irc_opt_status]
            if orig_status != new_status:
                if new_status == 'job_success':
                    update_field = {
                        irc_status: 'opt_success', irc_opt_status: new_status, 'irc_equal': 'waiting for check'
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)
                else:
                    update_field = {
                        irc_status: 'opt_fail', irc_opt_status: new_status
                    }
                    qm_collection.update_one(
                        target, {"$set": update_field}, True)
