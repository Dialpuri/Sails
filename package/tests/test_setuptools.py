import os

def test_setup_existance():
    assert os.path.exists("setup.py"), "Swap to setuptools before pushing!"