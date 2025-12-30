XRR Class
=========

The ``XRR`` class is a tool for processing X-ray Reflectivity data obtained with a 2D detector.
This class provides a robust pipeline for converting raw detector images into quantitative reflectivity curves.

Overview
--------

X-ray Reflectivity (XRR) measures the intensity of X-rays reflected from a surface as a function of the incidence angle.
The resulting reflectivity curve contains interference fringes (Kiessig fringes) that are characteristic of the Scattering Length Density (SLD) profile perpendicular to the surface.

The ``XRR`` class simplifies the analysis of ID10 data by handling the full reduction process:

*   **2D Data Reduction:** The class processes 2D detector images for each point in the XRR scan. It defines regions of interest (ROI) for the specular signal and the background, performing integration and subtraction to obtain the reflected intensity.
*   **Normalization:** Intensities are normalized by the direct beam monitor and transmission to account for flux variations and attenuators used during the scan.
*   **Corrections:**
    *   **Footprint Correction:** Corrects for the geometric effect where the beam footprint on the sample exceeds the sample size at low angles.
    *   **Attenuator Calibration:** Automatically corrects for discrepancies in transmission factors by analyzing overlapping data points (doubles) measured with different attenuator settings.
*   **Data Representation:** Calculates the scattering vector $q_z$ and generates the reflectivity curve $R(q_z)$.
*   **Visualization:** Offers methods to plot the reflectivity curve (log scale), :math:`R \cdot q_z^4`, and 2D reciprocal space maps (:math:`q_x` vs :math:`q_z`) to inspect off-specular scattering.
*   **ORSO Support:** Supports saving data in the ORSO (Open Reflectometry Standards Organization) format, ensuring long-term data preservation and interoperability.

Key Features
------------

*   **Automatic Attenuation Handling:** Detects and corrects for mismatches between attenuator steps, ensuring a smooth and continuous reflectivity curve.
*   **Footprint Correction:** Built-in geometric correction for finite sample sizes, crucial for accurate low-angle data.
*   **q-space Map Generation:** Visualizes the intensity distribution in ($q_x, q_z$) space, allowing for the detection of diffuse scattering or sample misalignment.
*   **Standardized Output:** Save results in simple ASCII formats or the rich ORSO standard format.

.. autoclass:: ESRF_ID10_SURF.XRR.XRR
    :members:
    :undoc-members:
    :show-inheritance:
