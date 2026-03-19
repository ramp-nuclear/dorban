from pathlib import Path

import pytest
from testbook import testbook


@pytest.mark.parametrize(
    "notebook",
    [
        "HTGR",
        "IAEA3D",
        "PWRMOX",
        "PWR_rod_bundle",
    ],
)
def test_notebook(notebook):
    d = (Path(__file__).parent / Path("../docs/source")).resolve()
    with testbook(str(d / f"{notebook}.ipynb"), execute=True) as tb:
        assert tb.ref("k")
