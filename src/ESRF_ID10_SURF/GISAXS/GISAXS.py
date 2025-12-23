import numpy as np


class GISAXS:
    def __init__(self, file, scans, detector_name='lambda',
                 alpha_i_name='chi', energy_name='monoe', sdd_name = 'sdd', **kwargs):
        """
        Initialize the GISAXS processor.
        Args:
            file (str): Path to the HDF5 file.
            scans (list or int): List of scan numbers or a single scan number.
            detector_name (str, optional): Name of detector dataset. Defaults to 'lambda'.
            alpha_i_name (str, optional): Name of incident angle dataset. Defaults to 'chi'.
            energy_name (str, optional): Name of the energy motor. Defaults to 'monoe'.
            sdd_name (str, optional): Name of the sample-to-detector distance motor. Defaults to 'sdd'.
            **kwargs: Arbitrary keyword arguments.
        """
        self.file = file
        self.scans = np.array(scans)
        self.detector_name = detector_name
        self.alpha_i_name = alpha_i_name
        self.energy_name = energy_name
        self.sdd_name = sdd_name
