import argparse
import sys
import os
import logging
import yaml

# Ensure we can import the package modules even if run as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    __package__ = "ESRF_ID10_SURF"

from .GID.GID import GID


def parse_scans(scan_str):
    """
    Parse a string of scans into a list of integers.

    Args:
        scan_str: String of scans to parse
    Returns:
        list of integers
    """
    scans = []
    parts = scan_str.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            scans.extend(range(start, end + 1))
        else:
            scans.append(int(part))
    return scans


def make_saving_dir(saving, filename):
    if saving == 'default':
        if 'RAW_DATA' in filename:
            file_dir, name_file = os.path.split(filename)
            saving_dir = file_dir.replace('RAW_DATA', 'SCRIPTS')
        else:
            logging.info(
                'Filename does not contain RAW_DATA. Files will be saved into the current working directory.')
            saving_dir = os.getcwd()
    else:
        saving_dir = saving
    return saving_dir


def main():
    parser = argparse.ArgumentParser(description="ESRF ID10 SURF Mythen 2 Calibration CLI")
    parser.add_argument("filename", help="Filename of scans to parse")
    parser.add_argument("scan_number", help="Scan number")
    args = parser.parse_args()

    if not os.path.exists(args.filename):
        print(f"Error: Data file '{args.config}' not found.")
        sys.exit(1)

    saving_dir = make_saving_dir("default", args.filename)
    saving_file = saving_dir+f"/Mythen2_calib_{args.scan_number}.png"

    try:
        os.makedirs(saving_dir, exist_ok=True)
    except OSError as e:
        print('Saving directory is impossible: ', e)

    PPD, mythen_gap = GID.calibrate_mythen(filename=args.filename, scanN=args.scan_number, plot=saving_file)

    print(f"Pixel per degree: {PPD}")
    print(f"Mythen gap: {mythen_gap}")

if __name__ == "__main__":
    main()
