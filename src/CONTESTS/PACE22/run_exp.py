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

RUN_TESTS = False

inst_dir = 'vcred-test-graphs'
output_root_dir = 'results'

solver_name = 'VCReducer'



# this program will run [thread_cnt] processes, each running TestsRunner, which runs tests_runner_threads processes,
# each of which calls the solver process (and CPSAT might use many workers...)
thread_cnt = 1 if not RUN_TESTS else 1
tests_runner_threads = 4 if not RUN_TESTS else 1

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
cpsat_threads = 8


def setNumThreads(t):
    global tests_runner_threads
    tests_runner_threads = t

def createCommands():
    global inst_dir, output_root_dir

    for alg in algorithms:
        cmd = getDefaultCommand()
        cmd += ' --run_name=vc-reducer'
        solver_params = ' --alg=' + alg + \
                        ' --time=120' + \
                        ' --rep=5' + \
                        ' --gran=1'
        cmd += ' --solver_params=\'' + solver_params + '\''
        all_tests_commands.append(cmd)

def createTestsCommands():



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
    print(f'{RUN_TESTS=} {thread_cnt=} {tests_runner_threads=} {all_input_files=}')

    try:
        if not os.path.exists(output_root_dir):
            os.makedirs(output_root_dir)
    except Exception as e:
        print(f"Error creating directory: {e}")


    if RUN_TESTS:
        solver_name = 'cpsat_exp_1_no_run'

    createTestsCommands()
    print(all_tests_commands)
    exit(1)

    if RUN_TESTS:
        print('#CAUTION! Taking only a fraction of all tests, just to test if it works as intended...')
        all_tests_commands = all_tests_commands[0::5]

    print(f"All {len(all_tests_commands)} commands to run:", *all_tests_commands, sep='\n\n', end='\n\n')

    p = multiprocessing.Pool(thread_cnt)
    dss = p.map(runTestForCommand, all_tests_commands, chunksize=1)
    createTablesAndRankings()