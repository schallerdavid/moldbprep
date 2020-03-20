import multiprocessing
import os
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import shutil
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
    counter = 0
    with open(file_path, 'r') as sdf_file:
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
            if '>  <' in line:
                properties.append(line.strip()[4:-1])
            elif '> <' in line:
                properties.append(line.strip()[3:-1])
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
    sdf_text += '\n$$$$\n'
    return sdf_text


def sdf_text_worker(merged_results, vendors, num_mols, start_time, mol_counter, fragment_counter, drug_like_counter,
                    big_counter, parent_fragment_collector, parent_drug_like_collector, parent_big_collector):
    fragment_collector, drug_like_collector, big_collector = [], [], []
    for index, row in merged_results.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        mol = Chem.AddHs(mol)
        mol.SetProp('_Name', row['smiles'])
        properties = {vendor: row[vendor] for vendor in vendors}
        molecular_weight = ExactMolWt(mol)
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


def write_sdf(merged_results, mols_per_file, output_path, vendors, num_processes):
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
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
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
                                                   big_counter, fragment_collector, drug_like_collector, big_collector))
                     for job_id in range(job_start, job_end)]
        for process in processes:
            process.start()
        for process in processes:
            process.join()
        if fragment_counter.value > mols_per_file:
            fragment.write('\n'.join(fragment_collector[0:mols_per_file - fragment_counter.value]))
            fragment_collector = fragment_collector[mols_per_file - fragment_counter.value:]
            fragment_counter.value = len(fragment_collector)
            fragment.close()
            fragment_file_counter += 1
            fragment = open(os.path.join(output_path, 'fragment_{}.sdf'.format(fragment_file_counter)))
        if len(fragment_collector) >= 1000:
            fragment.write('\n'.join(fragment_collector))
            fragment_collector = manager.list()
        if drug_like_counter.value > mols_per_file:
            drug_like.write('\n'.join(drug_like_collector[0:mols_per_file - drug_like_counter.value]))
            drug_like_collector = drug_like_collector[mols_per_file - drug_like_counter.value:]
            drug_like_counter.value = len(drug_like_collector)
            drug_like.close()
            drug_like_file_counter += 1
            drug_like = open(os.path.join(output_path, 'drug_like_{}.sdf'.format(drug_like_file_counter)))
        if len(drug_like_collector) >= 1000:
            drug_like.write('\n'.join(drug_like_collector))
            drug_like_collector = manager.list()
        if big_counter.value > mols_per_file:
            big.write('\n'.join(big_collector[0:mols_per_file - big_counter.value]))
            big_collector = big_collector[mols_per_file - big_counter.value:]
            big_counter.value = len(big_collector)
            big.close()
            big_file_counter += 1
            big = open(os.path.join(output_path, 'big_{}.sdf'.format(big_file_counter)))
        if len(big_collector) >= 1000:
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
