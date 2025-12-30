.. esrf_id10_surf documentation master file, created by
   sphinx-quickstart on Sun Dec 28 21:56:32 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: https://img.shields.io/pypi/v/esrf_id10_surf.svg
    :target: https://pypi.org/project/esrf_id10_surf/
    :alt: PyPI version

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT
    :alt: License: MIT

.. image:: https://img.shields.io/pypi/pyversions/esrf_id10_surf.svg
    :target: https://pypi.org/project/esrf_id10_surf/
    :alt: Python Versions

=================
ID10-SURF data processing
=================

Welcome to the documentation for ``esrf_id10_surf``, a Python package for analyzing surface X-ray scattering data from the ID10 beamline at the European Synchrotron Radiation Facility (ESRF).

This package provides tools for processing:

*   **Grazing Incidence Diffraction (GID)**: For analyzing in-plane structure.
*   **X-ray Reflectivity (XRR)**: For analyzing out-of-plane structure, thickness, and roughness.
*   **Grazing Incidence Small-Angle X-ray Scattering (GISAXS)**: (Under development)

------------------------
Installation
------------------------

To install the library create a virtual environment and run:

   pip install esrf-id10-surf

If you would like to install the library and modify the code to see the effects immediately:

   git clone https://github.com/chelberserker/ESRF_ID10_SURF
   cd ESRF_ID10_SURF
   pip install -e .

------------------------
Documentation
------------------------
.. toctree::
   :maxdepth: 2
   :caption: Documentation:

   GID
   XRR
   GISAXS

Indices and tables

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


------------------------
Usage examples
------------------------

This code is accompanied by a simple example of usage:

   :ref:`Examples`
