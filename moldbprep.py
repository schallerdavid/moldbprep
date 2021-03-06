"""
moldbprep.py
Prepare, standardize and merge molecule databases for virtual screening.

Handles the primary functions
"""
import argparse
import math
from moldbprep.io import count_sdf_mols, database_prompt, write_sdf, write_statistics, time_to_text
from moldbprep.standardize import standardize_mols, merge_ids
import multiprocessing
import os
import pandas as pd
import time

logo = '\n'.join(["                               _     _ _              v. alpha    ",
                  "               _ __ ___   ___ | | __| | |__  _ __  _ __ ___ _ __  ",
                  "              | '_ ` _ \ / _ \| |/ _` | '_ \| '_ \| '__/ _ \ '_ \ ",
                  "              | | | | | | (_) | | (_| | |_) | |_) | | |  __/ |_) |",
                  "              |_| |_| |_|\___/|_|\__,_|_.__/| .__/|_|  \___| .__/ ",
                  "                                            |_|            |_|    ",
                  "                Prepare, standardize and merge molecule databases ",
                  "                             for virtual screening.               "])


def standardize_processes(sdf_file_dict, mols_per_job):
    jobs = []
    for sdf_path, value in sdf_file_dict.items():
        num_mols = value[0]
        vendor = value[1]
        identifier_field = value[2]
        for mol_start in range(0, num_mols, mols_per_job):
            if mol_start + mols_per_job > num_mols:
                mol_end = num_mols - 1
            else:
                mol_end = mol_start + mols_per_job - 1
            jobs.append({'sdf_path': sdf_path, 'mol_start': mol_start, 'mol_end': mol_end, 'vendor': vendor,
                         'identifier_field': identifier_field})
    return jobs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='moldbprep', description=logo, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', dest='input_paths',
                        help='directory containing sdf files, or paths to input sdf files seperated by comma',
                        required=True)
    parser.add_argument('-o', dest='output_path', help='path to output folder', default='.')
    parser.add_argument('-p', dest='num_processes', help='number of parallel processes', default=1)
    parser.add_argument('-m', dest='mols_per_file', help='number of molecules per file', default=1000000)
    parser.add_argument('-s', dest='max_stereo_isomers',
                        help='maximal number of stereo isomers to generate per molecule, default=0 --> no enumeration',
                        default=0)
    parser.add_argument('-t', dest='tautomer', action='store_true', help='standardize tautomers')
    parser.add_argument('-e', dest='embed', action='store_true', help='generate 3D conformation')
    parser.add_argument('-a', dest='addhs', action='store_true', help='add hydrogens')
    parser.add_argument('-v', dest='verbose', action='store_true', help='show RDKit warnings')
    if os.path.isdir(parser.parse_args().input_paths):
        input_directory = os.path.abspath(parser.parse_args().input_paths)
        input_paths = [os.path.join(input_directory, path) for path in os.listdir(input_directory)
                       if path[-4:] == '.sdf']
    else:
        input_paths = [os.path.abspath(path) for path in parser.parse_args().input_paths.split(',')]
    output_path = os.path.abspath(parser.parse_args().output_path)
    num_processes = int(parser.parse_args().num_processes)
    mols_per_file = int(parser.parse_args().mols_per_file)
    max_stereo_isomers = int(parser.parse_args().max_stereo_isomers)
    tautomer = parser.parse_args().tautomer
    embed = parser.parse_args().embed
    addhs = parser.parse_args().addhs
    verbose = parser.parse_args().verbose
    print(logo)
    sdf_file_dict = {file_path: [count_sdf_mols(file_path), *database_prompt(file_path)] for file_path in input_paths}
    vendors = list(set([value[1] for value in sdf_file_dict.values()]))
    num_mols = sum([value[0] for value in sdf_file_dict.values()])
    max_mols_per_db = max([value[0] for value in sdf_file_dict.values()])
    print('Standardizing {} molecules from {} vendors...'.format(num_mols, len(vendors)))
    start_time = time.time()
    manager = multiprocessing.Manager()
    results = manager.list()
    failures = manager.list()
    jobs = manager.list()
    for job in standardize_processes(sdf_file_dict, math.ceil(max_mols_per_db / num_processes)):
        jobs.append(job)
    mol_counter = multiprocessing.Value('i', 0)
    processes = [multiprocessing.Process(target=standardize_mols, args=(jobs, mol_counter, num_mols, results,
                                                                        start_time, vendors, max_stereo_isomers,
                                                                        failures, tautomer, verbose))
                 for process_id in range(num_processes)]
    for process in processes:
        process.start()
    for process in processes:
        process.join()
    print('Processing results...')
    results = pd.DataFrame(list(results), columns=['smiles'] + vendors)
    results = merge_ids(results, vendors)
    print('Writing {} molecules...'.format(results.shape[0]))
    write_sdf(results, mols_per_file, output_path, vendors, failures, num_processes, addhs, embed, verbose)
    failures = list(failures)
    if len(failures) > 0:
        with open(os.path.join(output_path, 'moldbprep.failures'), 'w') as file:
            file.write('\n'.join(failures))
    write_statistics(num_mols, results, vendors, output_path, len(failures))
    print('Finished after {}.'.format(time_to_text(time.time() - start_time)))
