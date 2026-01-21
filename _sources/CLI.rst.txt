Command Line Interface
=======================

Command Line Interface (CLI) provides a convenient wrapper for batch processing of the XRR and GID data.
This interface reads a configuration YAML file with a following structure:

Instrument definition:

.. code-block:: yaml
    instrument: 'id10-surf'

Next, visit definition for traceability and metadata storage in ORSO format:

.. code-block:: yaml
    visit:
      local_contact: "John Doe"
      user: 'Big Prof'
      user_affiliation: 'Harvard'
      visit id: 'sc7777'
      date: 2025-12-26
      saving: 'default'     # can be either /path/to/the/dir or cwd. default is PROCESSED_DATA if path contains RAW_DATA, cwd otherwise


Next, we define an XRR setup by providing reduction parameters:

.. code-block:: yaml
    setup_xrr:
      PX0: 401
      PY0: 300
      dPX: 8
      dPY: 5
      bckg_gap: 3
      monitor_name: 'ionch2'
      alpha_i_name: 'mu'
      beam_size: 20
      sample_size: 17

Similar for GID reduction parameters:

.. code-block:: yaml
    setup_gid:
      PX0: 50
      PPD: 198.5
      mythen_gap: 90
      alpha_i_name: 'mu'
      monitor_name: 'ionch2'
      I0: 2e12

Next we define files and scans for processing XRR. It is possible to define multiple files

.. code-block:: yaml
    xrr:
      - file: '/data/visitor/ls3582/id10-surf/20251120/RAW_DATA/sample/sample_dataset/sample_dataset.h5'
        scans:
          - {zgH: 7, refl: "8,9"}
          - {zgH: 14, refl: "15,16"}
      - file: '/data/visitor/ls3582/id10-surf/20251120/RAW_DATA/lipids/lipids_40C/lipids_40C.h5'
        scans:
          - {zgH: 5, refl = "6,7"}

Similar definition is for GID processing. Note, that several scans on one line will stitch them together and those should not overlap.

.. code-block:: yaml
    gid:
      - file: '/data/ls3582/id10-surf/20251120/RAW_DATA/lipids/lipids_20C/lipids_20C.h5'
        scans:
          - 11
          - 19, 20


To call CLI simply run:

.. code-block:: bash
    python3 cli.py processing_config.yaml