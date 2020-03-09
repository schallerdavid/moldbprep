

def count_sdf_mols(file_path):
    """
    This function returns the number of molecules in an sdf-file.

    Parameters
    ----------
    file_path :  str
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
