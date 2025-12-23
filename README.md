# ESRF ID10 SURF Data Analysis

This package provides a collection of tools for analyzing surface X-ray scattering data from the ID10 beamline at the European Synchrotron Radiation Facility (ESRF).

## Features

*   **XRR (X-ray Reflectivity):** Process and analyze X-ray reflectivity data.
*   **GID (Grazing Incidence Diffraction):** Process and analyze grazing incidence diffraction data.
*   **GISAXS (Grazing Incidence Small-Angle X-ray Scattering):** (Under development)

## Installation

To install the package, clone the repository and install the dependencies:

```bash
git clone https://github.com/esrf-id10-surf/esrf-id10-surf.git
cd esrf-id10-surf
pip install -e .
```

## Usage

```python
from esrf_id10_surf.xrr import XRR

# Example usage
xrr_data = XRR(
    file='data.h5',
    scans=[1, 2, 3],
    # ... other parameters
)

xrr_data.plot_reflectivity()
```
