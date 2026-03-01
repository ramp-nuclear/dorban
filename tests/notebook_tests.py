import os
from pathlib import Path

import pytest
from testbook import testbook


def test_notebook_HTGR():
    p=Path.cwd()
    d = Path("../docs/source")
    os.chdir(d)
    with testbook(f"HTGR.ipynb",execute=True) as tb:
        assert tb.ref("k")
    os.chdir(p)


def test_notebook_IAEA3D():
    p=Path.cwd()
    d = Path("../docs/source")
    os.chdir(d)
    with testbook(f"IAEA3D.ipynb", execute=True) as tb:
        assert tb.ref("k")
    os.chdir(p)

@pytest.mark.slow
def test_notebook_PWRMOX():
    p=Path.cwd()
    d = Path("../docs/source")
    os.chdir(d)
    with testbook(f"PWRMOX.ipynb",execute=True) as tb:
        assert tb.ref("k")
    os.chdir(p)

def test_notebook_PWR_rod_bundle():
    p=Path.cwd()
    d = Path("../docs/source")
    os.chdir(d)
    with testbook(f"PWR_rod_bundle.ipynb", execute=True) as tb:
        assert tb.ref("k")
    os.chdir(p)
