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

inst_dir = 'testing_inputs' if RUN_TESTS else 'input_minimal'
output_root_dir = 'testing_results' if RUN_TESTS else 'results_minimal'

input_minimal_dir, results_minimal_dir = 'input_minimal', 'results_minimal'
input_large_dir, results_large_dir = 'input_large', 'results_large'
solver_name = 'cpsat_exp_1'



# this program will run [thread_cnt] processes, each running TestsRunner, which runs tests_runner_threads processes,
# each of which calls the solver process (and CPSAT might use many workers...)
thread_cnt = 1 if not RUN_TESTS else 1
tests_runner_threads = 1 if not RUN_TESTS else 1
single_core_solver_threads = 8

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

    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir
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


algorithms = ["ihs", "mtz", "diverses", "dreyfvs"]
cycle_enumeration_types = [1,2,3]
mtz_cycle_augmentation = [0,1,2]
pi_arcs_percentage = np.arange(0.05, 0.31, 0.05)
init_single_iter_time = [1,2,3]
cycle_scales = ['linear_h', 'linear', 'log', 'sqrt']
cpsat_threads = 8

large_time_sec = 900
small_time_sec = 300

# --threads=6
# --mtd=ex_mtd_file.csv
# --time=300
# --iter_time=1
# --find_optimal=false
# --alg=diverses
# --log_cpsat_progress=false
# --cycle_enumeration=1
# --init_L=4
# --max_new_cycles_iter_scale=linear
# --mtz_auxiliary_cycles_mode=2
# --pi_arcs_perc_to_add=0.1
# --use_init_sol_as_hint_mode=3

def setNumThreads(t):
    global tests_runner_threads
    tests_runner_threads = t

def createCycleEnumerationCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir

    for alg in algorithms:
        for cet in cycle_enumeration_types:
            cmd = getDefaultCommand()
            cmd += ' --run_name=cycle_enumeration__' + alg + '_' + str(cet)
            solver_params = ' --time=' + str(small_time_sec) + \
                            ' --threads=' + str(cpsat_threads) + \
                            ' --alg=' + alg + \
                            ' --cycle_enumeration=' + str(cet)
            cmd += ' --solver_params=\'' + solver_params + '\''
            all_tests_commands.append(cmd)

def createMTZCycleAugmentationCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir

    for ca in mtz_cycle_augmentation:
        cmd = getDefaultCommand()
        cmd += ' --run_name=mtz_cycle_augmentation__mtz_' + str(ca)
        solver_params = ' --time=' + str(small_time_sec) + \
                        ' --threads=' + str(cpsat_threads) + \
                        ' --alg=mtz' + \
                        ' --mtz_auxiliary_cycles_mode=' + str(ca)
        cmd += ' --solver_params=\'' + solver_params + '\''
        all_tests_commands.append(cmd)

def createCycleScalingCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir

    for alg in algorithms:
        setNumThreads( single_core_solver_threads if alg in ['diverses', 'dreyfvs'] else 1)

        for cs in cycle_scales:
            cmd = getDefaultCommand()
            cmd += ' --run_name=cycle_scaling__' + alg + '_' + cs
            solver_params = ' --time=' + str(small_time_sec) + \
                            ' --threads=' + str(cpsat_threads) + \
                            ' --alg=' + alg + \
                            ' --max_new_cycles_iter_scale=' + cs
            cmd += ' --solver_params=\'' + solver_params + '\''
            all_tests_commands.append(cmd)
    setNumThreads(1)

def createPiArcsPercentageCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir

    for alg in algorithms:
        setNumThreads( single_core_solver_threads if alg in ['diverses', 'dreyfvs'] else 1 )

        for pap in pi_arcs_percentage:
            cmd = getDefaultCommand()
            cmd += ' --run_name=pi_arcs_percentage__' + alg + '_' + str(pap)
            solver_params = ' --time=' + str(small_time_sec) + \
                            ' --threads=' + str(cpsat_threads) + \
                            ' --alg=' + alg + \
                            ' --pi_arcs_perc_to_add=' + str(pap)
            cmd += ' --solver_params=\'' + solver_params + '\''
            all_tests_commands.append(cmd)

    setNumThreads(1)

def createInitialSingleIterationTimeCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_minimal_dir, results_minimal_dir

    for alg in filter(lambda s : s.lower() != 'hs', algorithms):
        setNumThreads( single_core_solver_threads if alg in ['diverses', 'dreyfvs'] else 1)

        for sit in init_single_iter_time:
            cmd = getDefaultCommand()
            cmd += ' --run_name=init_single_iter_time__' + alg + '_' + str(sit)
            solver_params = ' --time=' + str(small_time_sec) + \
                            ' --threads=' + str(cpsat_threads) + \
                            ' --alg=' + alg + \
                            ' --iter_time=' + str(sit)
            cmd += ' --solver_params=\'' + solver_params + '\''
            all_tests_commands.append(cmd)

    setNumThreads(1)

def createAllGraphCommands():
    global inst_dir, output_root_dir
    inst_dir, output_root_dir = input_large_dir, results_large_dir

    for alg in algorithms:
        setNumThreads( single_core_solver_threads if alg in ['diverses', 'dreyfvs'] else 1)

        cmd = getDefaultCommand()
        cmd += ' --run_name=all_graphs__' + alg
        solver_params = ' --time=' + str(large_time_sec) + \
                        ' --threads=' + str(cpsat_threads) + \
                        ' --alg=' + alg
        cmd += ' --solver_params=\'' + solver_params + '\''
        all_tests_commands.append(cmd)

    setNumThreads(1)

def createTestsCommands():

    createCycleEnumerationCommands()
    createMTZCycleAugmentationCommands()
    createCycleScalingCommands()
    createPiArcsPercentageCommands()
    createInitialSingleIterationTimeCommands()
    createAllGraphCommands()


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
        if not os.path.exists(results_minimal_dir):
            os.makedirs(results_minimal_dir)
    except Exception as e:
        print(f"Error creating directory: {e}")

    try:
        if not os.path.exists(results_large_dir):
            os.makedirs(results_large_dir)
    except Exception as e:
        print(f"Error creating directory: {e}")

    try:
        if not os.path.exists(output_root_dir):
            os.makedirs(output_root_dir)
    except Exception as e:
        print(f"Error creating directory: {e}")


    if RUN_TESTS:
        solver_name = 'cpsat_exp_1_no_run'

    createTestsCommands()

    if RUN_TESTS:
        print('#CAUTION! Taking only a fraction of all tests, just to test if it works as intended...')
        all_tests_commands = all_tests_commands[0::5]

    print(f"All {len(all_tests_commands)} commands to run:", *all_tests_commands, sep='\n\n', end='\n\n')

    p = multiprocessing.Pool(thread_cnt)
    dss = p.map(runTestForCommand, all_tests_commands, chunksize=1)
    createTablesAndRankings()