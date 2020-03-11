"""
moldbprep.py
Prepare, standardize and merge molecule databases for virtual screening.

Handles the primary functions
"""
import argparse
import time

logo = '\n'.join(["                               _     _ _              v. alpha    ",
                  "               _ __ ___   ___ | | __| | |__  _ __  _ __ ___ _ __  ",
                  "              | '_ ` _ \ / _ \| |/ _` | '_ \| '_ \| '__/ _ \ '_ \ ",
                  "              | | | | | | (_) | | (_| | |_) | |_) | | |  __/ |_) |",
                  "              |_| |_| |_|\___/|_|\__,_|_.__/| .__/|_|  \___| .__/ ",
                  "                                            |_|            |_|    ",
                  "                Prepare, standardize and merge molecule databases ",
                  "                             for virtual screening.               "])

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(prog='moldbprep', description=logo, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', dest='input_paths', help='paths to input sdf files seperated by comma', required=True)
    parser.add_argument('-o', dest='output_path', help='path to output sdf file', required=True)
    input_paths = parser.parse_args().input_paths.split(',')
    output_path = parser.parse_args().output_path
    print(logo)
    print('Paths to input sdf files:\n', '\n '.join(input_paths))
    print('Path to output sdf file:\n', output_path)
    print('Finished after {} s.'.format(time.time() - start_time))
