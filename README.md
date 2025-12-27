# ESRF ID10-SURF Data Analysis

This package provides a collection of tools for analyzing surface X-ray scattering data from the ID10-SURF beamline at the European Synchrotron Radiation Facility (ESRF).

## Features

*   **XRR (X-ray Reflectivity):** Process and analyze X-ray reflectivity data obtained with a 2D pixel detector
*   **GID (Grazing Incidence Diffraction):** Process and analyze grazing incidence diffraction data taken with a linear detector
*   **GIS(W)AXS (Grazing Incidence Small(Wide)-Angle X-ray Scattering):** (Under development)

## Installation

To install the package use pip:

```bash
pip install esrf-id10-surf
```

## Usage
Detailed example is given in the jupyter notebook. 
```python
from ESRF_ID10_SURF.XRR import XRR

# Example usage
xrr_data = XRR(
    file='data.h5',
    scans=[1, 2, 3],
    # ... other parameters
)

xrr_data.apply_auto_corrections(sample_size=2, beam_size=10)
xrr_data.plot_reflectivity()
```
