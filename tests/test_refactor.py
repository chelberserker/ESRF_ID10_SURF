import pytest
import h5py
import numpy as np
import os
import logging
from ESRF_ID10_SURF.GID.GID import GID
from ESRF_ID10_SURF.XRR.XRR import XRR
from ESRF_ID10_SURF.GIXS.GIXS import GIXS

@pytest.fixture
def mock_h5_file(tmp_path):
    filepath = tmp_path / "test_data.h5"
    with h5py.File(filepath, 'w') as f:
        # Scan 1
        g1 = f.create_group("1.1")
        meas = g1.create_group("measurement")
        inst = g1.create_group("instrument/positioners")
        sample = g1.create_group("sample")

        # Datasets
        # GID/XRR/GIXS common
        meas.create_dataset("mon", data=np.ones(10))
        meas.create_dataset("autof_eh1_transm", data=np.ones(10)*0.1)
        meas.create_dataset("autof_eh1_curratt", data=np.zeros(10))
        meas.create_dataset("sec", data=np.ones(10))

        # Energy (usually in positioners)
        inst.create_dataset("monoe", data=12.4)

        # Sample name
        sample.create_dataset("name", data=b"TestSample") # Weird format emulation

        # Pi
        meas.create_dataset("fb_Pi", data=30.0)

        # GID specific
        # Detector (Mythen) 1D but GID loads it as... wait GID loads 2D data?
        # GID loads 'mythen2'. In __load_single_scan__, data = f.get(...)[()].
        # If mythen2 is 2D (points x pixels), usually.
        # GID code: nx, ny = np.shape(self.data).
        meas.create_dataset("mythen2", data=np.random.rand(10, 1280)) # 10 points, 1280 pixels
        meas.create_dataset("delta", data=np.linspace(0, 1, 10)) # Angle
        inst.create_dataset("chi", data=0.5) # alpha_i fixed for GID

        # XRR specific
        # Detector (mpx_cdte_22_eh1). Usually 3D (points x Y x X)?
        # XRR code: nic, nxc, nyc = np.shape(self.data).
        meas.create_dataset("mpx_cdte_22_eh1", data=np.random.rand(10, 200, 200))
        # alpha_i (chi or mu) scanning
        meas.create_dataset("chi_scan", data=np.linspace(0, 2, 10))

        # GIXS specific
        # Detector (eiger4m). 3D
        meas.create_dataset("eiger4m", data=np.random.rand(10, 100, 100))

        # Scan 2 (Append test)
        g2 = f.create_group("2.1")
        meas2 = g2.create_group("measurement")
        inst2 = g2.create_group("instrument/positioners")
        sample2 = g2.create_group("sample")

        meas2.create_dataset("mon", data=np.ones(5))
        meas2.create_dataset("autof_eh1_transm", data=np.ones(5)*0.1)
        meas2.create_dataset("autof_eh1_curratt", data=np.zeros(5))
        meas2.create_dataset("sec", data=np.ones(5))
        inst2.create_dataset("monoe", data=12.4)
        sample2.create_dataset("name", data=b"TestSample")
        meas2.create_dataset("fb_Pi", data=30.0)

        meas2.create_dataset("mythen2", data=np.random.rand(5, 1280))
        meas2.create_dataset("delta", data=np.linspace(1, 1.5, 5))
        inst2.create_dataset("chi", data=0.5)

        meas2.create_dataset("mpx_cdte_22_eh1", data=np.random.rand(5, 200, 200))
        meas2.create_dataset("chi_scan", data=np.linspace(2, 3, 5))

        meas2.create_dataset("eiger4m", data=np.random.rand(5, 100, 100))

    return str(filepath)

def test_gid_loading(mock_h5_file, tmp_path):
    # Test loading single scan
    gid = GID(file=mock_h5_file, scans=[1], saving_dir=str(tmp_path))
    assert gid.data.shape == (10, 1280)
    assert len(gid.monitor) == 10
    assert len(gid.angle) == 10
    assert gid.sample_name == "TestSample"

    # Test loading multiple scans
    gid_multi = GID(file=mock_h5_file, scans=[1, 2], saving_dir=str(tmp_path))
    # scan 1: 10 points. scan 2: 5 points. Skip points default=1.
    # Total = 10 + (5 - 1) = 14.
    assert gid_multi.data.shape[0] == 14
    assert len(gid_multi.monitor) == 14
    assert len(gid_multi.angle) == 14

def test_xrr_loading(mock_h5_file, tmp_path):
    # XRR uses alpha_i_name for scanning angle. In mock, 'chi_scan'.
    xrr = XRR(file=mock_h5_file, scans=[1], alpha_i_name='chi_scan', saving_dir=str(tmp_path))
    assert xrr.data.shape == (10, 200, 200)
    assert len(xrr.alpha_i) == 10

    xrr_multi = XRR(file=mock_h5_file, scans=[1, 2], alpha_i_name='chi_scan', saving_dir=str(tmp_path))
    assert xrr_multi.data.shape[0] == 14
    assert len(xrr_multi.alpha_i) == 14

def test_gixs_loading(mock_h5_file, tmp_path):
    # GIXS loading
    # GIXS expects a geometry object but it's optional in __init__?
    # GIXS init: geometry=None.
    # But integrate_single needs it. We only test loading here.
    gixs = GIXS(file=mock_h5_file, scan=[1], saving_dir=str(tmp_path))
    assert gixs.data.shape == (10, 100, 100)

    # Test multi scan (new feature)
    gixs_multi = GIXS(file=mock_h5_file, scan=[1, 2], saving_dir=str(tmp_path))
    assert gixs_multi.data.shape[0] == 14
