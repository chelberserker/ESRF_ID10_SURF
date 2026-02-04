
import pytest
import numpy as np
import h5py
from ESRF_ID10_SURF.XRR import XRR

def create_mock_data(filename, scan_num=1, n_points=10, img_size=516):
    with h5py.File(filename, "w") as f:
        scan_n = str(scan_num)
        base_path = f"{scan_n}.1"
        meas_path = f"{base_path}/measurement/"
        inst_path = f"{base_path}/instrument/positioners/"
        sample_path_str = f"{base_path}/sample/name"

        f.create_group(meas_path)
        f.create_group(inst_path)
        f.create_group(f"{base_path}/sample")

        rng = np.random.default_rng(42)
        # Note: In XRR.py, produce_Qmap uses range(516) for pixels, so we assume typical detector width.
        # But produce_Qmap hardcodes pixels = np.array(range(516)).
        # It doesn't use self.data shape for that part directly, but self.PX0.

        data = rng.integers(0, 100, size=(n_points, 100, 516)).astype(np.float64)
        # Note: __process_2D_data__ expects data to be shaped (N, Y, X).
        # We need X dimension to be compatible if we were using real data, but produce_Qmap logic relies on hardcoded 516.

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

        energy = np.array([12.4])
        f.create_dataset(f"{inst_path}monoe", data=energy)
        f.create_dataset(sample_path_str, data=b"TestSample")

@pytest.fixture
def mock_h5_file(tmp_path):
    filename = tmp_path / "test_qmap_data.h5"
    create_mock_data(filename, scan_num=1, n_points=10)
    return str(filename)

def test_produce_qmap(mock_h5_file):
    xrr = XRR(file=mock_h5_file, scans=[1])

    # Run produce_Qmap
    xrr.produce_Qmap()

    assert xrr.Qz_map is not None
    assert xrr.Qx_map is not None

    # Check shape: (N_points, 516)
    # n_points=10 in fixture
    assert xrr.Qz_map.shape == (10, 516)
    assert xrr.Qx_map.shape == (10, 516)

    # Check values are finite
    assert np.all(np.isfinite(xrr.Qz_map))
    assert np.all(np.isfinite(xrr.Qx_map))
