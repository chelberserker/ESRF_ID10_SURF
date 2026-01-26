import os
import logging
import time
import numpy as np
import matplotlib.pyplot as plt

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

class BaseSurf:
    """
    Base class for processing Surface X-ray Scattering data.
    """

    def __init__(self, file, scans,
                 alpha_i_name='chi', detector_name=None, monitor_name='mon',
                 transmission_name='autof_eh1_transm', att_name='autof_eh1_curratt',
                 cnttime_name='sec', energy_name='monoe',
                 I0=1e13, saving_dir=None, sample_name=None,
                 **kwargs):

        self.file = file
        self.scans = np.array(scans, ndmin=1)

        self.alpha_i_name = alpha_i_name
        self.detector_name = detector_name
        self.monitor_name = monitor_name
        self.transmission_name = transmission_name
        self.att_name = att_name
        self.cnttime_name = cnttime_name
        self.energy_name = energy_name

        self.I0 = I0
        self.saving_dir = saving_dir
        self.sample_name = sample_name or ""
        self.Pi = 100 # Default pressure

        # Data attributes
        self.data = np.empty((0, 0))
        self.alpha_i = None
        self.monitor = None
        self.transmission = None
        self.attenuator = None
        self.cnttime = None
        self.energy = None

        # List of attributes to append when loading multiple scans
        # Subclasses should extend this list
        self.appendable_keys = ['monitor', 'transmission', 'attenuator', 'cnttime']

        # Initialize other kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

    def _load_single_scan(self, scan_n):
        """
        Load data for a single scan.
        Must be implemented by subclasses.
        Returns:
            dict: Dictionary containing data arrays.
        """
        raise NotImplementedError("Subclasses must implement _load_single_scan")

    def load_data(self, skip_points=1):
        """
        Load data from all scans and concatenate.
        """
        t0 = time.time()
        logger.info("Start loading data.")

        # Load first scan
        first_scan_n = str(self.scans[0])
        first_scan_data = self._load_single_scan(first_scan_n)

        # Populate attributes from first scan
        for key, value in first_scan_data.items():
            setattr(self, key, value)

        # Load subsequent scans
        if len(self.scans) > 1:
            for scan_num in self.scans[1:]:
                scan_n = str(scan_num)
                try:
                    logger.info(f'Loading scan {scan_n}')
                    scan_data = self._load_single_scan(scan_n)

                    for key in self.appendable_keys:
                        if hasattr(self, key) and key in scan_data:
                            curr_val = getattr(self, key)
                            new_val = scan_data[key]

                            # Determine if we skip points
                            if new_val is not None:
                                try:
                                    # Check if new_val is a scalar (0-d array or float/int)
                                    if np.ndim(new_val) == 0:
                                        # It's a scalar (e.g. energy per scan), append without slicing
                                        if curr_val is None or np.size(curr_val) == 0:
                                            # If first scan was also scalar, curr_val is scalar.
                                            # We need to make it an array to append.
                                            # Note: if curr_val is None (unlikely if appendable), just set to new.
                                            # If curr_val is scalar, make array.
                                            if curr_val is not None:
                                                setattr(self, key, np.array([curr_val, new_val]))
                                            else:
                                                setattr(self, key, np.array([new_val]))
                                        else:
                                            # If curr_val is already array (from >1 scans), append
                                            setattr(self, key, np.append(curr_val, new_val))
                                    else:
                                        # It's an array (measurement points), apply skip_points
                                        sliced_val = new_val[skip_points:]
                                        if curr_val is None or np.size(curr_val) == 0:
                                            setattr(self, key, sliced_val)
                                        else:
                                            setattr(self, key, np.append(curr_val, sliced_val, axis=0))
                                except Exception as e:
                                    logger.warning(f"Could not append {key} for scan {scan_n}: {e}")

                    logger.info(f'Loaded scan #{scan_n}')
                except Exception as e:
                    logger.error(f"Error appending scan {scan_n}: {e}")

        logger.info("Loading completed. Reading time %3.3f sec" % (time.time() - t0))

    def _check_saving_dir(self):
        """
        Check if the saving directory is set; if not, set a default based on the current working directory and sample name.
        """
        if self.saving_dir:
            pass
        else:
            self.saving_dir = os.path.join(os.getcwd(), self.sample_name)

    def _ensure_sample_dir(self):
        """
        Ensure the saving directory exists, creating it if necessary.
        """
        if not self.saving_dir:
            self._check_saving_dir()

        try:
            # Basic validation: ensure path is valid
            os.makedirs(self.saving_dir, exist_ok=True)
        except OSError as e:
            logger.error('Saving directory is impossible: %s', e)

    def _save_figure(self, fig, suffix, dpi=200):
        """
        Helper method to save a matplotlib figure.
        """
        self._ensure_sample_dir()
        filename = os.path.join(self.saving_dir, f'{self.__class__.__name__}_{self.sample_name}_scan_{self.scans}_{suffix}.png')
        fig.savefig(filename, dpi=dpi)
        logger.info(f'Plot saved to {filename}')
