"""
moldbprep.py
Prepare, standardize and merge molecule databases for virtual screening.

Handles the primary functions
"""
import argparse
from moldbprep.io import count_sdf_mols, database_prompt
from moldbprep.standardize import standardize_mols
import multiprocessing
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


def generate_processes(sdf_file_dict, mols_per_job):
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
    parser.add_argument('-i', dest='input_paths', help='paths to input sdf files seperated by comma', required=True)
    parser.add_argument('-o', dest='output_path', help='path to output sdf file', required=True)
    parser.add_argument('-p', dest='num_processes', help='number of parallel processes', default=1)
    input_paths = parser.parse_args().input_paths.split(',')
    output_path = parser.parse_args().output_path
    num_processes = int(parser.parse_args().num_processes)
    print(logo)
    sdf_file_dict = {file_path: [count_sdf_mols(file_path), *database_prompt(file_path)] for file_path in input_paths}
    vendors = [value[1] for value in sdf_file_dict.values()]
    num_mols = sum([value[0] for value in sdf_file_dict.values()])
    start_time = time.time()
    manager = multiprocessing.Manager()
    results = manager.list()
    jobs = manager.list()
    for job in generate_processes(sdf_file_dict, 500):
        jobs.append(job)
    mol_counter = multiprocessing.Value('i', 0)
    processes = [multiprocessing.Process(target=standardize_mols, args=(jobs, mol_counter, num_mols, results,
                                                                        start_time, vendors)) for process_id in
                 range(num_processes)]
    for process in processes:
        process.start()
    for process in processes:
        process.join()
    results = pd.DataFrame(list(results), columns=['smiles', 'inchikey'] + vendors)
    print('Finished after {} s.'.format(time.time() - start_time))
