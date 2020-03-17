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


def sdf_text(mol, mol_name, properties):
    mol.SetProp('_Name', mol_name)
    text = Chem.MolToMolBlock(mol)
    text += '\n'.join(['>  <{}>\n{}\n'.format(key, value) for key, value in properties.items()])
    text += '\n$$$$\n'
    return text


def write_sdf(merged_results, mols_per_file, output_path):
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
    mol_counter = 0
    num_mols = merged_results.shape[0]
    fragment = open(os.path.join(output_path, 'fragment.sdf'), 'w')
    fragment_counter = 0
    fragment_collector = []
    drug_like = open(os.path.join(output_path, 'drug-like.sdf'), 'w')
    drug_like_counter = 0
    drug_like_collector = []
    big = open(os.path.join(output_path, 'big.sdf'), 'w')
    big_counter = 0
    big_collector = []
    column_names = list(merged_results.columns.values)
    # os.mkdir(os.path.join(output_path, 'covalent'))
    # covalent_counter = 0
    # os.mkdir(os.path.join(output_path, 'building-block'))
    # building_block_counter_counter = 0
    start_time = time.time()
    for index, row in merged_results.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        mol = Chem.AddHs(mol)
        properties = {name: row[name] for name in column_names}
        molecular_weight = ExactMolWt(mol)
        if molecular_weight < 1200:
            if molecular_weight < 300:
                fragment_counter += 1
                fragment_collector.append(sdf_text(mol, 'fragment_' + str(fragment_counter), properties))
                if fragment_counter % mols_per_file == 0 or len(fragment_collector) >= 1000:
                    fragment.write('\n'.join(fragment_collector))
                    fragment_collector = []
                    if fragment_counter % mols_per_file == 0:
                        fragment.close()
                        fragment = open(os.path.join(output_path, 'fragment_{}.sdf'.format(
                            int(fragment_counter / mols_per_file) + 1)))
            elif 300 <= molecular_weight < 700:
                drug_like_counter += 1
                drug_like_collector.append(sdf_text(mol, 'drug_like_' + str(drug_like_counter), properties))
                if drug_like_counter % mols_per_file == 0 or len(drug_like_collector) >= 1000:
                    drug_like.write('\n'.join(drug_like_collector))
                    drug_like_collector = []
                    if drug_like_counter % mols_per_file == 0:
                        drug_like.close()
                        drug_like = open(os.path.join(output_path, 'drug_like_{}.sdf'.format(
                            int(drug_like_counter / mols_per_file) + 1)))
            else:
                big_counter += 1
                big_collector.append(sdf_text(mol, 'big_' + str(big_counter), properties))
                if big_counter % mols_per_file == 0 or len(big_collector) >= 1000:
                    big.write('\n'.join(big_collector))
                    big_collector = []
                    if big_counter % mols_per_file == 0:
                        big.close()
                        big = open(os.path.join(output_path, 'big_{}.sdf'.format(
                            int(big_counter / mols_per_file) + 1)))
        mol_counter += 1
        update_progress(mol_counter / num_mols, 'Progress of writing',
                        ((time.time() - start_time) / mol_counter) * (num_mols - mol_counter))
    if len(fragment_collector) > 0:
        fragment.write('\n'.join(fragment_collector))
        fragment.close()
    if len(drug_like_collector) > 0:
        drug_like.write('\n'.join(drug_like_collector))
        drug_like.close()
    if len(big_collector) > 0:
        big.write('\n'.join(big_collector))
        big.close()
    return
