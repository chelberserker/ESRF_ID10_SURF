
import h5py
import numpy as np
import pytest
import os
import shutil
from ESRF_ID10_SURF.XRR import XRR

def create_mock_data(filename, scan_num=1, n_points=50, img_size=516):
    with h5py.File(filename, "w") as f:
        scan_n = str(scan_num)
        base_path = f"{scan_n}.1"
        meas_path = f"{base_path}/measurement/"
        inst_path = f"{base_path}/instrument/positioners/"
        sample_path_str = f"{base_path}/sample/name"

        # Create groups
        f.create_group(meas_path)
        f.create_group(inst_path)
        # Ensure parent group exists for sample/name
        f.create_group(f"{base_path}/sample")

        # Create datasets
        # Data: (points, Y, X)
        # Using a fixed seed for reproducibility in tests
        rng = np.random.default_rng(42)
        data = rng.integers(0, 100, size=(n_points, img_size, img_size)).astype(np.float64)
        f.create_dataset(f"{meas_path}mpx_cdte_22_eh1", data=data)

        alpha_i = np.linspace(0, 1, n_points)
        f.create_dataset(f"{meas_path}chi", data=alpha_i)

        monitor = rng.random(n_points) * 1000
        f.create_dataset(f"{meas_path}mon", data=monitor)

        transmission = np.ones(n_points)
        f.create_dataset(f"{meas_path}autof_eh1_transm", data=transmission)

        attenuator = np.zeros(n_points)
        f.create_dataset(f"{meas_path}autof_eh1_curratt", data=attenuator)

        cnttime = np.ones(n_points)
        f.create_dataset(f"{meas_path}sec", data=cnttime)

        f.create_dataset(f"{meas_path}fb_Pi", data=np.array([100.0]))

        energy = np.array([12.4]) # 12.4 keV
        f.create_dataset(f"{inst_path}monoe", data=energy)

        # sample name
        f.create_dataset(sample_path_str, data=b"TestSample")

@pytest.fixture
def mock_h5_file(tmp_path):
    filename = tmp_path / "test_xrr_data.h5"
    create_mock_data(filename, scan_num=1, n_points=50, img_size=100) # Smaller size for tests
    return str(filename)

def test_xrr_processing(mock_h5_file):
    xrr = XRR(file=mock_h5_file, scans=[1])

    # Check attributes are populated
    assert xrr.reflectivity is not None
    assert xrr.reflectivity.shape == (50,)
    assert xrr.Smap2D is not None
    assert xrr.Smap2D.shape == (50, 100) # Smap2D has size N x X_pixels

    # Simple sanity checks
    assert np.all(np.isfinite(xrr.reflectivity))
    assert np.all(np.isfinite(xrr.Smap2D))

    # Check cleanup happens implicitly by temp dir removal, or we can explicit check nothing is left if we were creating side files.
    # XRR creates directories if saving data. Here we don't call save methods.
