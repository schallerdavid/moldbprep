

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
