import argparse
import sys
import os
import logging
import yaml

# Ensure we can import the package modules even if run as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    __package__ = "ESRF_ID10_SURF"

from .XRR.XRR import XRR
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


def main():
    parser = argparse.ArgumentParser(description="ESRF ID10 SURF Data Processing CLI")
    parser.add_argument("config", help="Path to the YAML configuration file")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Configuration file '{args.config}' not found.")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading YAML config: {e}")
        sys.exit(1)

    # Defaults from setup sections
    setup_xrr = config.get('setup_xrr', {})
    setup_gid = config.get('setup_gid', {})
    visit_conf = config.get('visit', {})

    user = visit_conf.get('user', 'ESRF')
    owner = visit_conf.get('user_affiliation', 'Your_University')
    saving = visit_conf.get('saving', 'default')

    
    # Check for XRR batch processing
    if 'xrr' in config and isinstance(config['xrr'], list):
        print("Starting XRR Batch Processing...")
        for item in config['xrr']:
            filename = item.get('file')
            if not filename:
                print("Warning: Skipped XRR item without 'file' key.")
                continue
            
            if not os.path.exists(filename):
                 print(f"Warning: File '{filename}' not found. Skipping.")
                 continue

            if saving == 'default':
                if 'RAW_DATA' in filename:
                    file_dir, name_file = os.path.split(filename)
                    saving_dir = file_dir.replace('RAW_DATA', 'PROCESSED_DATA')
                else:
                    logging.info(
                        'Filename does not contain RAW_DATA. Files will be saved into the current working directory.')
                    saving_dir = os.getcwd()
            else:
                saving_dir = saving

            scans_list = item.get('scans', [])
            for scan_item in scans_list:
                try:
                    # Parse parameters
                    alpha_i = setup_xrr.get('alpha_i_name', 'mu')
                    monitor = setup_xrr.get('monitor_name', 'ionch2')
                    px0 = int(setup_xrr.get('PX0', 404))
                    py0 = int(setup_xrr.get('PY0', 165))
                    dpx = int(setup_xrr.get('dPX', 3))
                    dpy = int(setup_xrr.get('dPY', 3))
                    bckg_gap = int(setup_xrr.get('bckg_gap', 1))
                    
                    z_scan = None
                    refl_scans = []

                    # determine scan numbers and normalization
                    if isinstance(scan_item, dict):
                         # Format: {zgH: 7, refl: "8,9"}
                         z_val = scan_item.get('zgH')
                         refl_val = scan_item.get('refl')
                         
                         if z_val:
                             z_scan = parse_scans(str(z_val)) # returns list
                         if refl_val:
                             refl_scans = parse_scans(str(refl_val))
                    else:
                        # Simple scan number or range string
                        refl_scans = parse_scans(str(scan_item))

                    print(f"Processing XRR params: file={filename}, scans={refl_scans}, z_scan={z_scan}")

                    if not refl_scans:
                        continue

                    processor = XRR(filename, refl_scans, alpha_i_name=alpha_i, monitor_name=monitor, PX0=px0, PY0=py0, dPX=dpx, dPY=dpy, bckg_gap=bckg_gap, saving_dir=saving_dir)

                    # Apply normalization if z-scan exists
                    if z_scan:
                        try:
                             z_processor = XRR(filename, z_scan, alpha_i_name='zgH', monitor_name=monitor, PX0=px0, PY0=py0)

                             sample_size = float(setup_xrr.get('sample_size', 1.0))
                             beam_size = float(setup_xrr.get('beam_size', 10.0))
                             
                             processor.apply_auto_corrections(sample_size=sample_size, beam_size=beam_size, z_scan=z_processor)
                        except Exception as e:
                            print(f"Error applying auto-corrections: {e}")


                    # For now, always save .dat
                    processor.save_reflectivity(format='dat')

                    processor.save_reflectivity(format='orso', owner=owner, creator=user)
                    processor.plot_reflectivity(save=True)
                    
                except Exception as e:
                    print(f"Error processing XRR scan {scan_item} in {filename}: {e}")

    # Check for GID batch processing
    if 'gid' in config and isinstance(config['gid'], list):
         print("Starting GID Batch Processing...")
         for item in config['gid']:
            filename = item.get('file')
            if not filename or not os.path.exists(filename):
                print(f"Warning: File '{filename}' not found or missing. Skipping.")
                continue


            if saving == 'default':
                if 'RAW_DATA' in filename:
                    file_dir, name_file = os.path.split(filename)
                    saving_dir = file_dir.replace('RAW_DATA', 'PROCESSED_DATA')
                else:
                    logging.info(
                        'Filename does not contain RAW_DATA. Files will be saved into the current working directory.')
                    saving_dir = os.getcwd()
            else:
                saving_dir = saving

            scans_list = item.get('scans', [])

            for scan_item in scans_list:
                try:
                    # Parse parameters
                    alpha_i = setup_gid.get('alpha_i_name', 'mu')
                    monitor = setup_gid.get('monitor_name', 'ionch2')
                    px0 = int(setup_gid.get('PX0', 404))
                    ppd = int(setup_gid.get('PPD', 165))
                    mythen_gap = int(setup_gid.get('mythen_gap', 100))
                    I0 = float(setup_gid.get('I0', 1e12))

                    gid_scans = []

                    # determine scan number
                    try:
                        gid_scans = parse_scans(str(scan_item))
                    except Exception as e:
                        print(f"Error parsing GID scan {scan_item}: {e}")

                    try:
                        gid_processor = GID(filename, gid_scans, PX0=px0, PPD=ppd, mythen_gap=mythen_gap, alpha_i_name=alpha_i, monitor_name=monitor, I0=I0, saving_dir=saving_dir)
                        gid_processor.plot_quick_analysis(save=True)
                        gid_processor.save_image_h5()
                    except Exception as e:
                        logging.error(f"Error processing GID scan {scan_item}: {e}")

                except Exception as e:
                    logging.error(f"Error parsing GID scan setup parameters: {e}")


if __name__ == "__main__":
    main()
