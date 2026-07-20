import csv
import json
import os
import sys
import shutil
import concurrent.futures
import multiprocessing
import timeit
import time
import platform
import argparse
import numpy as np
import itertools
from pathlib import Path

inst_dir = 'vcred-exp-graphs'
output_root_dir = 'results'

# inst_dir = 'vcred-tests-small'
# output_root_dir = 'results-tests-small'

solver_name = 'VCReducer'



# this program will run [thread_cnt] processes, each running TestsRunner, which runs tests_runner_threads processes,
# each of which calls the solver process (and CPSAT might use many workers...)
thread_cnt = 2
tests_runner_threads = 2

def getDefaultCommand():
    cmd = 'python3 TestsRunner.py' + \
          ' --instances_dir=' + inst_dir + \
          ' --output_root_dir=' + output_root_dir + \
          ' --threads=' + str(tests_runner_threads) + \
          ' --solver_name=' + solver_name + \
          ' --remove_existing_results=false' + \
          ' --skip_existing_results=true' + \
          ' --compute=true' + \
          ' --run_judge=false' + \
          ' --create_rankings=false' + \
          ' --rerun_failed_tests=true' + \
          ' --require_metadata_creation=true'
    return cmd

def createTablesAndRankings():
    global inst_dir, output_root_dir

    cmd = 'python3 TestsRunner.py' + \
          ' --instances_dir=' + inst_dir + \
          ' --output_root_dir=' + output_root_dir + \
          ' --compute=true' + \
          ' --run_solver=false' + \
          ' --run_judge=false' + \
          ' --create_verdict_file=false' + \
          ' --create_rankings=true' + \
          ' --threads=1' + \
          ' --solver_name=' + solver_name + \
          ' --remove_existing_results=false' + \
          ' --skip_existing_results=true' + \
          ' --report_runs_in_separate_lines=true'

    print('\nCreating tables and ranking, running command', cmd)
    os.system(cmd)

all_tests_commands = []


algorithms = ["cpsat-sat", "cpsat-def", "numvc"]


def setNumThreads(t):
    global tests_runner_threads
    tests_runner_threads = t

def createTestsCommands():
    global inst_dir, output_root_dir

    for def1 in [True,False]:
        for alg in algorithms:
            setNumThreads( (16 // thread_cnt) if alg == 'numvc' else (4 // thread_cnt) )
            solver_time = ( 60 if alg == 'numvc' else 120 )

            cmd = getDefaultCommand()
            cmd += ' --run_name=vc-reducer_def1-' + str(def1)
            solver_params = '--alg=' + alg + \
                            ' --time=' + str(solver_time) + \
                            ' --rep=5' + \
                            ' --gran=1' + \
                            ' --use_def1_dom=' + str(def1)
            cmd += ' --solver_params=\'' + solver_params + '\''
            all_tests_commands.append(cmd)


    print(f'\nThere are altogether {len(all_tests_commands)} commands to run in total')

def runTestForCommand(cmd):
    print('Running command', cmd)
    os.system(cmd)

def count_files(root_dir, extensions):
    root = Path(root_dir)
    extensions = {ext.lower() for ext in extensions}  # normalize

    return sum(
        1 for f in root.rglob("*")
        if f.is_file() and f.suffix.lower() in extensions
    )

# Example


if __name__ == '__main__':
    all_input_files = count_files(inst_dir, ['.txt', '.in', '.mtx', '.edges'])
    print(all_input_files)

    print(f'{(platform.system())=}')
    print(f'{thread_cnt=} {tests_runner_threads=} {all_input_files=}')

    try:
        if not os.path.exists(output_root_dir):
            os.makedirs(output_root_dir)
    except Exception as e:
        print(f"Error creating directory: {e}")


    createTestsCommands()

    print(f"All {len(all_tests_commands)} commands to run:", *all_tests_commands, sep='\n\n', end='\n\n')

    p = multiprocessing.Pool(thread_cnt)
    dss = p.map(runTestForCommand, all_tests_commands, chunksize=1)
    createTablesAndRankings()