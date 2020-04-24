import multiprocessing
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.Chem.Descriptors import ExactMolWt
import re
import sys
import time


def count_sdf_mols(file_path):
    """
    This function returns the number of molecules in an sdf-file.

    Parameters
    ----------
    file_path : str
        Full path to sdf file.

    Returns
    -------
    counter : int
        Number of molecules.

    """
    print('Counting molecules in {}...'.format(file_path))
    counter = 0
    with open(file_path, 'r', errors='backslashreplace') as sdf_file:
        for line in sdf_file:
            if '$$$$' in line:
                counter += 1
    return counter


def sdf_properties(file_path):
    """
    This function returns a list of properties stored in an sdf-file.

    Parameters
    ----------
    file_path : str
         Full path to sdf file.

    Returns
    -------
    properties : list
        Properties stored in sdf-file.

    """
    properties = []
    with open(file_path, 'r') as sdf_file:
        for line in sdf_file:
            if '>  <' in line or '> <' in line:
                properties.append(re.search('<.*>|$', line).group()[1:-1])
            elif '$$$$' in line:
                break
    return properties


def database_prompt(file_path):
    """
    This function prompts the user to enter the vendor name and to identify the sdf field storing the molecule
    identifier in an sdf file.

    Parameters
    ----------
    file_path : str
         Full path to sdf file.

    Returns
    -------
    vendor : str
        Name of vendor.

    identifier_field : str
        Name of sdf field storing the molecule identifier.

    """
    vendor = ''
    id_column = 0
    properties = ['None'] + sdf_properties(file_path)
    while len(vendor) < 1:
        vendor = input('Provide a vendor name for sdf file located at {}.\n>>> '.format(file_path))
    while id_column not in range(1, len(properties) + 1):
        id_column = int(input('Enter the number for the sdf field storing the molecule identifier.\n' + '\n'.join(
            '{} - '.format(counter + 1) + property for counter, property in enumerate(properties)) + '\n>>> '))
    identifier_field = properties[int(id_column) - 1]
    return vendor, identifier_field


def time_to_text(seconds):
    """
    This function converts a time in seconds into a reasonable format.

    Parameters
    ----------
    seconds : float
        Time in seconds.

    Returns
    -------
    time_as_text: str
        Time in s, min, h, d, weeks or years depending on input.

    """
    if seconds > 60:
        if seconds > 3600:
            if seconds > 86400:
                if seconds > 1209600:
                    if seconds > 62899252:
                        time_as_text = 'years'
                    else:
                        time_as_text = '{} weeks'.format(round(seconds / 1209600, 1))
                else:
                    time_as_text = '{} d'.format(round(seconds / 86400, 1))
            else:
                time_as_text = '{} h'.format(round(seconds / 3600, 1))
        else:
            time_as_text = '{} min'.format(round(seconds / 60, 1))
    else:
        time_as_text = '{} s'.format(int(seconds))
    return time_as_text


def update_progress(progress, progress_info, eta):
    """
    This function writes a progress bar to the terminal.

    Parameters
    ----------
    progress: float
        Progress of process described by number between 0 and 1.

    progress_info: str
        Info text that should be placed before the progress bar.

    eta: float
        Estimated time needed for finishing the process.

    """
    bar_length = 10
    block = int(bar_length * progress)
    if progress == 1.0:
        status = '         Done\n'
    else:
        status = '  ETA {:8}'.format(time_to_text(eta))
    text = '\r{}: [{}] {:>5.1f}%{}'.format(progress_info, '=' * block + ' ' * (bar_length - block), progress * 100,
                                           status)
    sys.stdout.write(text)
    sys.stdout.flush()
    return


def sdf_text(mol, properties):
    """
    This function converts an RDKit molecule into an sdf representation as text.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        An RDKit molecule.

    mol_name: str
        Name of the molecule.

    properties: dict
        Dictionary of sdf properties with property name as key and property value as value.

    Returns
    -------
    sdf_text: str
        Molecule as text in sdf format.

    """
    sdf_text = Chem.MolToMolBlock(mol)
    sdf_text += '\n'.join(['>  <{}>\n{}\n'.format(key, value) for key, value in properties.items()])
    sdf_text += '\n$$$$'
    return sdf_text


def sdf_text_worker(merged_results, vendors, num_mols, start_time, mol_counter, fragment_counter, drug_like_counter,
                    big_counter, parent_fragment_collector, parent_drug_like_collector, parent_big_collector,
                    failures, addhs, embed, verbose):
    if not verbose:
        RDLogger.DisableLog('rdApp.*')
    fragment_collector, drug_like_collector, big_collector = [], [], []
    for index, row in merged_results.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if addhs:
                mol = Chem.AddHs(mol)
            if embed:
                AllChem.EmbedMolecule(mol)
            mol.SetProp('_Name', row['smiles'])
            properties = {vendor: row[vendor] for vendor in vendors}
            molecular_weight = ExactMolWt(mol)
        except:
            failures.append(' '.join(['write_error', row['smiles']]))
            molecular_weight = 10000
        if molecular_weight < 1200:
            if molecular_weight < 300:
                with fragment_counter.get_lock():
                    fragment_counter.value += 1
                fragment_collector.append(sdf_text(mol, properties))
            elif 300 <= molecular_weight < 700:
                with drug_like_counter.get_lock():
                    drug_like_counter.value += 1
                drug_like_collector.append(sdf_text(mol, properties))
            else:
                with big_counter.get_lock():
                    big_counter.value += 1
                big_collector.append(sdf_text(mol, properties))
        with mol_counter.get_lock():
            mol_counter.value += 1
            update_progress(mol_counter.value / num_mols, 'Progress of writing',
                            ((time.time() - start_time) / mol_counter.value) * (num_mols - mol_counter.value))
    parent_fragment_collector.extend(fragment_collector)
    parent_drug_like_collector.extend(drug_like_collector)
    parent_big_collector.extend(big_collector)
    return


def write_sdf(merged_results, mols_per_file, output_path, vendors, failures, num_processes, addhs, embed, verbose):
    """
    This function writes molecules to sd-files with vendor IDs as properties.

    Parameters
    ----------
    merged_results: pandas.DataFrame
        Processed molecules with smiles and vendor IDs.

    mols_per_file: int
        Number of molecules writter per file.

    output_path: str
        Directory for writing sd-files.

    """
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    manager = multiprocessing.Manager()
    mol_counter = multiprocessing.Value('i', 0)
    num_mols = merged_results.shape[0]
    fragment = open(os.path.join(output_path, 'fragment.sdf'), 'w')
    fragment_file_counter = 0
    fragment_counter = multiprocessing.Value('i', 0)
    fragment_collector = manager.list()
    drug_like = open(os.path.join(output_path, 'drug-like.sdf'), 'w')
    drug_like_file_counter = 0
    drug_like_counter = multiprocessing.Value('i', 0)
    drug_like_collector = manager.list()
    big = open(os.path.join(output_path, 'big.sdf'), 'w')
    big_file_counter = 0
    big_counter = multiprocessing.Value('i', 0)
    big_collector = manager.list()
    mols_per_job = 1000
    mols_in_memory = 10000
    jobs = []
    for mol_start in range(0, num_mols, mols_per_job):
        if mol_start + mols_per_job <= num_mols:
            jobs.append((mol_start, mol_start + mols_per_job))
        else:
            jobs.append((mol_start, num_mols))
    job_chunks = []
    for job_start in range(0, len(jobs), num_processes):
        if job_start + num_processes <= len(jobs):
            job_chunks.append((job_start, job_start + num_processes))
        else:
            job_chunks.append((job_start, len(jobs)))
    start_time = time.time()
    for job_start, job_end in job_chunks:
        processes = [multiprocessing.Process(target=sdf_text_worker,
                                             args=(merged_results[jobs[job_id][0]: jobs[job_id][1]], vendors, num_mols,
                                                   start_time, mol_counter, fragment_counter, drug_like_counter,
                                                   big_counter, fragment_collector, drug_like_collector, big_collector,
                                                   failures, addhs, embed, verbose))
                     for job_id in range(job_start, job_end)]
        for process in processes:
            process.start()
        for process in processes:
            process.join()
        if fragment_counter.value > mols_per_file:
            fragment.write('\n'.join(fragment_collector[0:mols_per_file - fragment_counter.value]))
            fragment_collector = manager.list(fragment_collector[mols_per_file - fragment_counter.value:])
            fragment_counter.value = len(fragment_collector)
            fragment.close()
            fragment_file_counter += 1
            fragment = open(os.path.join(output_path, 'fragment_{}.sdf'.format(fragment_file_counter)), 'w')
        if len(fragment_collector) >= mols_in_memory:
            fragment.write('\n'.join(fragment_collector))
            fragment_collector = manager.list()
        if drug_like_counter.value > mols_per_file:
            drug_like.write('\n'.join(drug_like_collector[0:mols_per_file - drug_like_counter.value]))
            drug_like_collector = manager.list(drug_like_collector[mols_per_file - drug_like_counter.value:])
            drug_like_counter.value = len(drug_like_collector)
            drug_like.close()
            drug_like_file_counter += 1
            drug_like = open(os.path.join(output_path, 'drug-like_{}.sdf'.format(drug_like_file_counter)), 'w')
        if len(drug_like_collector) >= mols_in_memory:
            drug_like.write('\n'.join(drug_like_collector))
            drug_like_collector = manager.list()
        if big_counter.value > mols_per_file:
            big.write('\n'.join(big_collector[0:mols_per_file - big_counter.value]))
            big_collector = manager.list(big_collector[mols_per_file - big_counter.value:])
            big_counter.value = len(big_collector)
            big.close()
            big_file_counter += 1
            big = open(os.path.join(output_path, 'big_{}.sdf'.format(big_file_counter)), 'w')
        if len(big_collector) >= mols_in_memory:
            big.write('\n'.join(big_collector))
            big_collector = manager.list()
    if len(fragment_collector) > 0:
        fragment.write('\n'.join(fragment_collector))
    if len(drug_like_collector) > 0:
        drug_like.write('\n'.join(drug_like_collector))
    if len(big_collector) > 0:
        big.write('\n'.join(big_collector))
    fragment.close()
    drug_like.close()
    big.close()
    return


def write_statistics(num_mols, merged_results, vendors, output_path, failure_count):
    """
    Write statistics about merged databases.

    Parameters
    ----------
    num_mols : int
        Number of input molecules.
    merged_results : pandas.DataFrame
        Dataframe containing the merged results.

    vendors : list
        List of vendors.

    output_path : str
        Path to output directory.

    failure_count : int
        Number of failures during standardization and writing.

    """
    vendor_matches = {vendor: merged_results[vendor] != '' for vendor in vendors}
    with open(os.path.join(output_path, 'database.statistics'), 'w') as file:
        file.write('Input: {} molecules\n\n'.format(num_mols))
        file.write('Vendor\tTotal\tUnique\n')
        for vendor in vendors:
            total = vendor_matches[vendor].sum()
            if len(vendors) > 1:
                unique = (vendor_matches[vendor] > pd.concat([vendor_matches[x] for x in vendors if x != vendor], axis=1
                                                             ).max(axis=1)).sum()
            else:
                unique = total
            file.write('\t'.join([vendor, str(total), str(unique)]) + '\n')
        file.write('\nCategory\tTotal\n')
        directory_contents = os.listdir(output_path)
        for file_name in ['fragment', 'drug-like', 'big']:
            mol_count = 0
            for directory_content in directory_contents:
                if file_name in directory_content:
                    mol_count_file = count_sdf_mols(os.path.join(output_path, directory_content))
                    mol_count += mol_count_file
                    if mol_count_file == 0:
                        os.remove(os.path.join(output_path, file_name + '.sdf'))
            file.write('{}\t{}\n'.format(file_name, mol_count))
        file.write('\nfailures\t{}'.format(failure_count))
    return
