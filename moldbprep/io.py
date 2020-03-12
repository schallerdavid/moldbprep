import sys


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
    seconds : int
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
