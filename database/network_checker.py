from connect import db
import subprocess
import os
from os import path
import sys
from openbabel import pybel
from helper import highlight_text, print_header


def select_ard_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    reg_query = {"ard_status":
                 {"$in":
                  ["job_launched", "job_running", "job_queueing"]
                  }
                 }
    targets = list(qm_collection.find(reg_query))

    return targets


def check_ard_status(job_id):
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


def check_ard_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ard_target()

    qm_collection = db['qm_calculate_center']

    # 2. check the job pbs status
    for target in targets:
        job_id = target['ard_jobid']
        # 2. check the job pbs status
        new_status = check_ard_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = 'job_finished'

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ard_status']
        if orig_status != new_status:

            if new_status == 'job_finished':
                update_field = {
                    'ard_status': new_status
                }
            else:
                update_field = {
                    'ard_status': new_status
                }

            qm_collection.update_one(target, {"$set": update_field}, True)


def print_information(generations):
    """
    For a given generations (int) print the information in database
    """
    qm_collection = db['qm_calculate_center']
    reactions_collection = db['reactions']

    gen_query = {"generations":
                 {"$in":
                  [generations]
                  }
                 }
    ard_query_0 = {'$and':
                   [
                       {"ard_status":
                        {"$in":
                         ["job_unrun"]}
                        },
                       {'generations': generations - 1}
                   ]
                   }
    ard_query_1 = {'$and':
                   [
                       {"ard_status":
                        {"$in":
                         ["job_launched", "job_running", "job_queueing"]}
                        },
                       {'generations': generations - 1}
                   ]
                   }
    ard_query_2 = {'$and':
                   [
                       {"ard_status":
                        {"$in":
                         ["job_finished"]}
                        },
                       {'generations': generations}
                   ]
                   }
    ssm_query_0 = {'$and':
                   [
                       {"ssm_status":
                        {"$in":
                         ["job_unrun"]}
                        },
                       {'generations': generations}
                   ]
                   }
    ssm_query_1 = {'$and':
                   [
                       {"ssm_status":
                        {"$in":
                         ["job_running", "job_queueing", "job_launched"]}
                        },
                       {'generations': generations}
                   ]
                   }
    ssm_query_2 = {'$and':
                   [
                       {"ssm_status":
                        {"$in":
                         ["job_success"]}
                        },
                       {'generations': generations}
                   ]
                   }
    ssm_query_3 = {'$and':
                   [
                       {"ssm_status":
                        {"$in":
                         ['job_fail', 'total dissociation', 'Exiting early', 'Ran out of nodes', 'All uphill', 'ts_guess high energy']}
                        },
                       {'generations': generations}
                   ]
                   }
    ts_refine_query_0 = {'$and':
                         [
                             {"ts_refine_status":
                              {"$in":
                               ["job_unrun"]}
                              },
                             {'generations': generations}
                         ]
                         }
    ts_refine_query_1 = {'$and':
                         [
                             {"ts_refine_status":
                              {"$in":
                               ["job_running", "job_queueing", "job_launched"]}
                              },
                             {'generations': generations}
                         ]
                         }
    ts_refine_query_2 = {'$and':
                         [
                             {"ts_refine_status":
                              {"$in":
                               ['job_success']}
                              },
                             {'generations': generations}
                         ]
                         }
    ts_refine_query_3 = {'$and':
                         [
                             {"ts_refine_status":
                              {"$in":
                               ['job_fail', 'All positive frequency', 'Have more than one imaginary frequency', 'Imaginary frequency greater than -50 cm-1']}
                              },
                             {'generations': generations}
                         ]
                         }
    ts_query_0 = {'$and':
                  [
                      {"ts_status":
                       {"$in":
                        ['job_unrun']}
                       },
                      {'generations': generations}
                  ]
                  }
    ts_query_1 = {'$and':
                  [
                      {"ts_status":
                       {"$in":
                        ['job_running', "job_queueing", "job_launched"]}
                       },
                      {'generations': generations}
                  ]
                  }
    ts_query_2 = {'$and':
                  [
                      {"ts_status":
                       {"$in":
                        ['job_success']}
                       },
                      {'generations': generations}
                  ]
                  }
    ts_query_3 = {'$and':
                  [
                      {"ts_status":
                       {"$in":
                        ['job_fail', 'All positive frequency', 'Have more than one imaginary frequency', 'Imaginary frequency greater than -50 cm-1']}
                       },
                      {'generations': generations}
                  ]
                  }
    irc_query_0 = {'$and':
                   [
                       {"irc_status":
                        {"$in":
                         ['job_unrun']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_1 = {'$and':
                   [
                       {"irc_status":
                        {"$in":
                         ['job_running', "job_queueing", "job_launched"]}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_2 = {'$and':
                   [
                       {"irc_status":
                        {"$in":
                         ['job_success']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_3 = {'$and':
                   [
                       {"irc_status":
                        {"$in":
                         ['job_fail']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_4 = {'$and':
                   [
                       {"irc_equal":
                        {"$in":
                         ['forward equal to reactant',
                          'backward equal to reactant']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_5 = {'$and':
                   [
                       {"irc_equal":
                        {"$in":
                         ['forward equal to reverse', 
                         'same forward and reverse reactant part but different active site',
                          'unintended']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_7 = {'$and':
                   [    {"irc_opt_status":'job_success'},
                       {"irc_equal":
                        {"$nin":
                         ['forward equal to reverse', 
                         'same forward and reverse reactant part but different active site',
                          'unintended',
                          'forward equal to reactant',
                          'backward equal to reactant']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_query_6 = {'$and':
                   [
                       {"irc_equal":
                        {"$in":
                         ['waiting for check']}
                        },
                       {'generations': generations}
                   ]
                   }
    irc_opt_query_0 = {'$and':
                       [
                           {"irc_opt_status":
                            {"$in":
                             ["job_unrun"]}
                            },
                           {'generations': generations}
                       ]
                       }
    irc_opt_query_1 = {'$and':
                       [
                           {"irc_forward_opt_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    irc_opt_query_4 = {'$and':
                       [
                           {"irc_backward_opt_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    irc_opt_query_2 = {'$and':
                       [
                           {"irc_opt_status":
                            {"$in":
                             ['job_success']}
                            },
                           {'generations': generations}
                       ]
                       }
    irc_opt_query_3 = {'$and':
                       [
                           {"irc_forward_opt_status":
                            {"$in":
                             ['job_fail', 'Have negative frequency']}
                            },
                           {'generations': generations}
                       ]
                       }
    irc_opt_query_5 = {'$and':
                       [
                           {"irc_backward_opt_status":
                            {"$in":
                             ['job_fail', 'Have negative frequency']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_0 = {'$and':
                       [
                           {"qmmm_opt_status":
                            {"$in":
                             ["job_unrun"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_1 = {'$and':
                       [
                           {"qmmm_opt_reactant_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_4 = {'$and':
                       [
                           {"qmmm_opt_product_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_2 = {'$and':
                       [
                           {"qmmm_opt_status":
                            {"$in":
                             ['job_success']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_3 = {'$and':
                       [
                           {"qmmm_opt_reactant_status":
                            {"$in":
                             ['job_fail', 'the qmmm ts fail, so delete this side']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_opt_query_5 = {'$and':
                       [
                           {"qmmm_opt_product_status":
                            {"$in":
                             ['job_fail', 'the qmmm ts fail, so delete this side']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_0 = {'$and':
                       [
                           {"qmmm_freq_status":
                            {"$in":
                             ["job_unrun"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_1 = {'$and':
                       [
                           {"qmmm_freq_reactant_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_4 = {'$and':
                       [
                           {"qmmm_freq_product_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_2 = {'$and':
                       [
                           {"qmmm_freq_status":
                            {"$in":
                             ['job_success']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_3 = {'$and':
                       [
                           {"qmmm_freq_reactant_status":
                            {"$in":
                             ['job_fail', 'Have negative frequency']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_query_5 = {'$and':
                       [
                           {"qmmm_freq_product_status":
                            {"$in":
                             ['job_fail', 'Have negative frequency']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_ts_query_0 = {'$and':
                       [
                           {"qmmm_freq_ts_status":
                            {"$in":
                             ["job_unrun", "restart"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_ts_query_1 = {'$and':
                       [
                           {"qmmm_freq_ts_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_ts_query_2 = {'$and':
                       [
                           {"qmmm_freq_ts_status":
                            {"$in":
                             ["job_success"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_ts_query_3 = {'$and':
                       [
                           {"qmmm_freq_ts_status":
                            {"$in":
                             ["job_fail"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_ts_freq_query_0 = {'$and':
                       [
                           {"qmmm_ts_freq_status":
                            {"$in":
                             ["job_unrun"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_ts_freq_query_1 = {'$and':
                       [
                           {"qmmm_ts_freq_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_ts_freq_query_2 = {'$and':
                       [
                           {"qmmm_ts_freq_status":
                            {"$in":
                             ["job_success"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_ts_freq_query_3 = {'$and':
                       [
                           {"qmmm_ts_freq_status":
                            {"$in":
                             ["job_fail", "Have more than one imaginary frequency", "All positive frequency"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_0 = {'$and':
                       [
                           {"qmmm_freq_opt_status":
                            {"$in":
                             ["job_unrun"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_1 = {'$and':
                       [
                           {"qmmm_freq_opt_reactant_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_4 = {'$and':
                       [
                           {"qmmm_freq_opt_product_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_2 = {'$and':
                       [
                           {"qmmm_freq_opt_status":
                            {"$in":
                             ['job_success']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_3 = {'$and':
                       [
                           {"qmmm_freq_opt_reactant_status":
                            {"$in":
                             ['job_fail', 'the qmmm ts fail, so delete this side']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_5 = {'$and':
                       [
                           {"qmmm_freq_opt_product_status":
                            {"$in":
                             ['job_fail', 'the qmmm ts fail, so delete this side']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_6 = {'$and':
                       [
                           {"qmmm_freq_opt_reactant_status":
                            {"$in":
                             ["restart"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_7 = {'$and':
                       [
                           {"qmmm_freq_opt_product_status":
                            {"$in":
                             ["restart"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_8 = {'$and':
                       [
                           {"qmmm_freq_opt_reactant_status":
                            {"$in":
                             ['qmmm ts fail']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_freq_opt_query_9 = {'$and':
                       [
                           {"qmmm_freq_opt_product_status":
                            {"$in":
                             ['qmmm ts fail']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_0 = {'$and':
                        [{"qmmm_sp_status":
                            {"$in":
                                ["job_unrun"]}
                            },
                            {"qmmm_sp_ts_status":
                            {"$in":
                                ["job_unrun"]}
                            },
                            {"generations":
                            {"$in":
                                [generations]}
                            }]}
    qmmm_sp_query_1 = {'$and':
                       [
                           {"qmmm_sp_reactant_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_4 = {'$and':
                       [
                           {"qmmm_sp_product_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_2 = {'$and':
                       [
                           {"qmmm_sp_status":
                            {"$in":
                             ['job_success']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_3 = {'$and':
                       [
                           {"qmmm_sp_reactant_status":
                            {"$in":
                             ['job_fail']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_5 = {'$and':
                       [
                           {"qmmm_sp_product_status":
                            {"$in":
                             ['job_fail']}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_6 = {'$and':
                       [
                           {"qmmm_sp_ts_status":
                            {"$in":
                             ["job_running", "job_queueing", "job_launched"]}
                            },
                           {'generations': generations}
                       ]
                       }
    qmmm_sp_query_7 = {'$and':
                       [
                           {"qmmm_sp_ts_status":
                            {"$in":
                             ['job_fail']}
                            },
                           {'generations': generations}
                       ]
                       }

    gen_targets = list(qm_collection.find(gen_query))
    ard_targets_0 = list(qm_collection.find(ard_query_0))
    ard_targets_1 = list(qm_collection.find(ard_query_1))
    ard_targets_2 = list(qm_collection.find(ard_query_2))
    ssm_targets_0 = list(qm_collection.find(ssm_query_0))
    ssm_targets_1 = list(qm_collection.find(ssm_query_1))
    ssm_targets_2 = list(qm_collection.find(ssm_query_2))
    ssm_targets_3 = list(qm_collection.find(ssm_query_3))
    ts_refine_targets_0 = list(qm_collection.find(ts_refine_query_0))
    ts_refine_targets_1 = list(qm_collection.find(ts_refine_query_1))
    ts_refine_targets_2 = list(qm_collection.find(ts_refine_query_2))
    ts_refine_targets_3 = list(qm_collection.find(ts_refine_query_3))
    ts_targets_0 = list(qm_collection.find(ts_query_0))
    ts_targets_1 = list(qm_collection.find(ts_query_1))
    ts_targets_2 = list(qm_collection.find(ts_query_2))
    ts_targets_3 = list(qm_collection.find(ts_query_3))
    irc_targets_0 = list(qm_collection.find(irc_query_0))
    irc_targets_1 = list(qm_collection.find(irc_query_1))
    irc_targets_2 = list(qm_collection.find(irc_query_2))
    irc_targets_3 = list(qm_collection.find(irc_query_3))
    irc_targets_4 = list(qm_collection.find(irc_query_4))
    irc_targets_5 = list(qm_collection.find(irc_query_5))
    irc_targets_6 = list(qm_collection.find(irc_query_6))
    irc_opt_targets_0 = list(qm_collection.find(irc_opt_query_0))
    irc_opt_targets_1 = list(qm_collection.find({'$or':
                                                    [irc_opt_query_1, irc_opt_query_4]
                                                    }))
    irc_opt_targets_2 = list(qm_collection.find(irc_opt_query_2))
    irc_opt_targets_3 = list(qm_collection.find({'$or':
                                                    [irc_opt_query_3, irc_opt_query_5]
                                                    }))
    invalid_product_target = list(qm_collection.find(irc_query_7))
    qmmm_opt_targets_0 = list(qm_collection.find(qmmm_opt_query_0))
    qmmm_opt_targets_1 = list(qm_collection.find({'$or':
                                                    [qmmm_opt_query_1, qmmm_opt_query_4]
                                                    }))
    qmmm_opt_targets_2 = list(qm_collection.find(qmmm_opt_query_2))
    qmmm_opt_targets_3 = list(qm_collection.find({'$or':
                                                    [qmmm_opt_query_3, qmmm_opt_query_5]
                                                    }))
    qmmm_freq_targets_0 = list(qm_collection.find(qmmm_freq_query_0))
    qmmm_freq_targets_1 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_query_1, qmmm_freq_query_4]
                                                    }))
    qmmm_freq_targets_2 = list(qm_collection.find(qmmm_freq_query_2))
    qmmm_freq_targets_3 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_query_3, qmmm_freq_query_5]
                                                    }))
    qmmm_freq_ts_targets_0 = list(qm_collection.find(qmmm_freq_ts_query_0))
    qmmm_freq_ts_targets_1 = list(qm_collection.find(qmmm_freq_ts_query_1))
    qmmm_freq_ts_targets_2 = list(qm_collection.find(qmmm_freq_ts_query_2))
    qmmm_freq_ts_targets_3 = list(qm_collection.find(qmmm_freq_ts_query_3))
    qmmm_ts_freq_targets_0 = list(qm_collection.find(qmmm_ts_freq_query_0))
    qmmm_ts_freq_targets_1 = list(qm_collection.find(qmmm_ts_freq_query_1))
    qmmm_ts_freq_targets_2 = list(qm_collection.find(qmmm_ts_freq_query_2))
    qmmm_ts_freq_targets_3 = list(qm_collection.find(qmmm_ts_freq_query_3))
    qmmm_freq_opt_targets_0 = list(qm_collection.find(qmmm_freq_opt_query_0))
    qmmm_freq_opt_targets_1 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_opt_query_1, qmmm_freq_opt_query_4]
                                                    }))
    qmmm_freq_opt_targets_2 = list(qm_collection.find(qmmm_freq_opt_query_2))
    qmmm_freq_opt_targets_3 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_opt_query_3, qmmm_freq_opt_query_5]
                                                    }))
    qmmm_freq_opt_targets_6 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_opt_query_6, qmmm_freq_opt_query_7]
                                                    }))
    qmmm_freq_opt_targets_8 = list(qm_collection.find({'$or':
                                                    [qmmm_freq_opt_query_8, qmmm_freq_opt_query_9]
                                                    }))
    qmmm_sp_targets_0 = list(qm_collection.find(qmmm_sp_query_0))
    qmmm_sp_targets_1 = list(qm_collection.find({'$or':
                                                    [qmmm_sp_query_1, qmmm_sp_query_4, qmmm_sp_query_6]
                                                    }))
    qmmm_sp_targets_2 = list(qm_collection.find(qmmm_sp_query_2))
    qmmm_sp_targets_3 = list(qm_collection.find({'$or':
                                                    [qmmm_sp_query_3, qmmm_sp_query_5, qmmm_sp_query_7]
                                                    }))
    reactant_target = list(reactions_collection.find({'generations': generations}))
    smi = []
    for target in reactant_target:
        reactant_smi = target['reactant_smiles']
        smi.append(reactant_smi)
    smi = set(smi)
    print(highlight_text('Generations : {}'.format(generations)))
    print('Reactant SMILES:')
    for smiles in smi:
        print(smiles)
    print('Nodes: {}'.format(len(gen_targets)))
    print(highlight_text('Unrun jobs'))
    if len(ard_targets_0) != 0:
        print('\n{} nodes should run ARD'.format(len(ard_targets_0)))
    if len(ssm_targets_0) != 0:
        print('{} nodes should run SSM'.format(len(ssm_targets_0)))
    if len(ts_refine_targets_0) != 0:
        print('{} nodes should run TS REFINE'.format(len(ts_refine_targets_0)))
    if len(ts_targets_0) != 0:
        print('{} nodes should run TS'.format(len(ts_targets_0)))
    if len(irc_targets_0) != 0:
        print('{} nodes should run IRC'.format(len(irc_targets_0)))
    if len(irc_opt_targets_0) != 0:
        print('{} nodes should run IRC OPT job'.format(len(irc_opt_targets_0)))
    if len(irc_targets_6) != 0:
        print('{} nodes are waiting for checking IRC EQUAL'.format(len(irc_targets_6)))
    if len(qmmm_opt_targets_0) != 0:
        print('{} nodes should run QMMM OPT job'.format(len(qmmm_opt_targets_0)))
    if len(qmmm_freq_opt_targets_0) != 0:
        print('{} nodes should run QMMM FREQ OPT job'.format(len(qmmm_freq_opt_targets_0)))
    if len(qmmm_freq_opt_targets_6) != 0:
        print('{} nodes should restart QMMM FREQ OPT job'.format(len(qmmm_freq_opt_targets_6)))
    if len(qmmm_freq_ts_targets_0) != 0:
        print('{} nodes should run QMMM TS job'.format(len(qmmm_freq_ts_targets_0)))
    if len(qmmm_freq_targets_0) != 0:
        print('{} nodes should run QMMM FREQ job'.format(len(qmmm_freq_targets_0)))
    if len(qmmm_ts_freq_targets_0) != 0:
        print('{} nodes should run QMMM TS FERQ job'.format(len(qmmm_ts_freq_targets_0)))
    if len(qmmm_sp_targets_0) != 0:
        print('{} nodes should run QMMM SP job\n'.format(len(qmmm_sp_targets_0)))
    print(highlight_text('Running jobs'))
    if len(ard_targets_1) != 0:
        print('\n{} nodes are running or queueing ARD'.format(len(ard_targets_1)))
    if len(ssm_targets_1) != 0:
        print('{} nodes are running or queueing SSM'.format(len(ssm_targets_1)))
    if len(ts_refine_targets_1) != 0:
        print('{} nodes are running or queueing TS REFINE'.format(len(ts_refine_targets_1)))
    if len(ts_targets_1) != 0:
        print('{} nodes are running or queueing TS'.format(len(ts_targets_1)))
    if len(irc_targets_1) != 0:
        print('{} nodes are running or queueing IRC'.format(len(irc_targets_1)))
    if len(irc_opt_targets_1) != 0:
        print('{} nodes are running or queueing in IRC OPT job'.format(len(irc_opt_targets_1)))
    if len(qmmm_opt_targets_1) != 0:
        print('{} nodes are running or queueing in QMMM OPT job'.format(len(qmmm_opt_targets_1)))
    if len(qmmm_freq_opt_targets_1) != 0:
        print('{} nodes are running, queueing in QMMM FREQ OPT job'.format(len(qmmm_freq_opt_targets_1)))
    if len(qmmm_freq_ts_targets_1) != 0:
        print('{} nodes are running or queueing QMMM TS'.format(len(qmmm_freq_ts_targets_1)))
    if len(qmmm_freq_targets_1) != 0:
        print('{} nodes are running or queueing QMMM FREQ'.format(len(qmmm_freq_targets_1)))
    if len(qmmm_ts_freq_targets_1) != 0:
        print('{} nodes are running or queueing QMMM TS FREQ'.format(len(qmmm_ts_freq_targets_1)))
    if len(qmmm_sp_targets_1) != 0:
        print('{} nodes are running or queueing in QMMM SP job\n'.format(len(qmmm_sp_targets_1)))
    print(highlight_text('Success jobs'))
    if len(ard_targets_2) != 0:
        print('\n{} nodes are success in ARD'.format(len(ard_targets_2)))
    if len(ssm_targets_2) != 0:
        print('{} nodes are success in SSM'.format(len(ssm_targets_2)))
    if len(ts_refine_targets_2) != 0:
        print('{} nodes are success in TS REFINE'.format(len(ts_refine_targets_2)))
    if len(ts_targets_2) != 0:
        print('{} nodes are success in TS'.format(len(ts_targets_2)))
    if len(irc_targets_2) != 0:
        print('{} nodes are success in IRC'.format(len(irc_targets_2)))
    if len(irc_opt_targets_2) != 0:
        print('{} nodes are success in IRC OPT job'.format(len(irc_opt_targets_2)))
    if len(qmmm_opt_targets_2) != 0:
        print('{} nodes are success in QMMM OPT job'.format(len(qmmm_opt_targets_2)))
    if len(qmmm_freq_opt_targets_2) != 0:
        print('{} nodes are success in QMMM FREQ OPT job'.format(len(qmmm_freq_opt_targets_2)))
    if len(qmmm_freq_ts_targets_2) != 0:
        print('{} nodes are success in QMMM TS job'.format(len(qmmm_freq_ts_targets_2)))
    if len(qmmm_freq_targets_2) != 0:
        print('{} nodes are success in QMMM FREQ job'.format(len(qmmm_freq_targets_2)))
    if len(qmmm_ts_freq_targets_2) != 0:
        print('{} nodes are success in QMMM TS FREQ job'.format(len(qmmm_ts_freq_targets_2)))
    if len(qmmm_sp_targets_2) != 0:
        print('{} nodes are success in QMMM SP job\n'.format(len(qmmm_sp_targets_2)))
    print(highlight_text('Failed jobs'))
    if len(ssm_targets_3) != 0:
        print('\n{} nodes are failed in SSM'.format(len(ssm_targets_3)))
    if len(ts_refine_targets_3) != 0:
        print('{} nodes are failed in TS REFINE'.format(len(ts_refine_targets_3)))
    if len(ts_targets_3) != 0:
        print('{} nodes are failed in TS'.format(len(ts_targets_3)))
    if len(irc_targets_3) != 0:    
        print('{} nodes are failed in IRC'.format(len(irc_targets_3)))
    if len(irc_opt_targets_3) != 0:
        print('{} nodes are failed in IRC OPT job'.format(len(irc_opt_targets_3)))
    if len(qmmm_opt_targets_3) != 0:
        print('{} nodes are failed in QMMM OPT job'.format(len(qmmm_opt_targets_3)))
    if len(qmmm_freq_opt_targets_3) != 0:
        print('{} nodes are failed in QMMM FREQ OPT job'.format(len(qmmm_freq_opt_targets_3)))
    if len(qmmm_freq_opt_targets_8) != 0:
        print('{} nodes are failed in one side'.format(len(qmmm_freq_opt_targets_8)))
    if len(qmmm_freq_ts_targets_3) != 0:
        print('{} nodes are failed in QMMM TS'.format(len(qmmm_freq_ts_targets_3)))
    if len(qmmm_freq_targets_3) != 0:
        print('{} nodes are failed in QMMM FREQ'.format(len(qmmm_freq_targets_3)))
    if len(qmmm_ts_freq_targets_3) != 0:
        print('{} nodes are failed in QMMM TS FREQ'.format(len(qmmm_ts_freq_targets_3)))
    if len(qmmm_sp_targets_3) != 0:
        print('{} nodes are failed in QMMM SP job\n'.format(len(qmmm_sp_targets_3)))
    print(highlight_text('IRC intended products'))
    if len(irc_targets_4) != 0:
        print('\n{} nodes are intended'.format(len(irc_targets_4)))
    if len(irc_targets_5) != 0:
        print('{} nodes are unintended'.format(len(irc_targets_5)))
    if len(invalid_product_target) != 0:
        print('{} products are invalid\n'.format(len(invalid_product_target)))

def update_network_status():
    status_collection = db['status']
    qm_collection = db['qm_calculate_center']
    statistics_collection = db['statistics']
    config_collection = db['config']
    target = list(config_collection.find({'generations': 1}))[0]
    qmmm = target['use_qmmm']
    use_irc = target['use_irc']
    if qmmm == '1':
        qmmm = True
    else:
        qmmm = False
    if use_irc == '1':
        use_irc = True
    else:
        use_irc = False

    ard_had_add_number = qm_collection.count_documents({})
    ard_should_add_number = sum([i['add how many products'] for i in list(statistics_collection.find({}))])

    ard_query = {"ard_status":
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
                         {'irc_backward_opt_status':
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

    qmmm_query = {'$or':
                     [
                         {"qmmm_opt_status":
                          {"$in":
                           ["job_unrun"]}},
                         {'qmmm_opt_reactant_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_opt_product_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_freq_opt_status':
                             {'$in':
                              ["job_unrun"]}},
                         {'qmmm_freq_opt_reactant_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_freq_opt_product_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_freq_status':
                             {'$in':
                              ["job_unrun"]}},
                         {'qmmm_freq_reactant_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_freq_product_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_freq_ts_status':
                             {'$in':
                              ["job_unrun", "job_launched", "job_running", "job_queueing"]}},
                         {'$and':[{'qmmm_sp_status':
                             {'$in':
                              ["job_unrun"]}},
                         {'qmmm_sp_ts_status':
                             {'$in':
                              ["job_unrun"]}}]},
                         {'qmmm_sp_reactant_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_sp_product_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}},
                         {'qmmm_sp_ts_status':
                             {'$in':
                              ["job_launched", "job_running", "job_queueing"]}}
                     ]
                     }

    if not use_irc:
        not_finished_number = len(list(qm_collection.find({'$or':
                                                [ssm_query, ts_refine_query, ts_query, insert_reaction_query, ard_query]
                                                })))
    else:
        if qmmm:
            not_finished_number = len(list(qm_collection.find({'$or':
                                                    [ssm_query, ts_refine_query, ts_query, irc_query, irc_opt_query, irc_equal_query, insert_reaction_query, ard_query, qmmm_query]
                                                    })))
        else:
            not_finished_number = len(list(qm_collection.find({'$or':
                                                    [ssm_query, ts_refine_query, ts_query, irc_query, irc_opt_query, irc_equal_query, insert_reaction_query, ard_query]
                                                    })))
             
    if ard_had_add_number == ard_should_add_number and not_finished_number == 0:
        print('Network converged')

        target = list(status_collection.find({}))
        update_field = {'status': 'Network converged'}
        status_collection.update_one(target[0], {"$set": update_field}, True)

print_header()
check_ard_jobs()

qm_collection = db['qm_calculate_center']
max_gen = qm_collection.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']
for i in range(max_gen):
    print_information(i+1)

update_network_status()
