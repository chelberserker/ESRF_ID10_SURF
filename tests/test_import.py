import unittest

def test_import_gid():
    from ESRF_ID10_SURF.GID import GID

def test_import_gixs():
    from ESRF_ID10_SURF.GIXS import GIXS

def test_import_xrr():
    from ESRF_ID10_SURF.XRR import XRR


try:
    test_import_gixs()
    test_import_xrr()
    test_import_gid()
except ImportError as e:
    print(e)
