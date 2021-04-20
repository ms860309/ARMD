from connect import db
import os
import fnmatch

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def get_time(file):
    with open(file) as f:
        log = f.read().splitlines()
    for line in reversed(log):
        if 'It took' in line:
            return int(line.split()[-2])

qm_collection = db['qm_calculate_center']
a = list(qm_collection.find({'generations':1}))

for i in a:
    dir_path = i['path']
    try:
        ssm_path = os.path.join(dir_path, 'SSM')
        ssm_output = find('ssm.job.o*', ssm_path)
        ssm_time = get_time(ssm_output)
    except:
        continue