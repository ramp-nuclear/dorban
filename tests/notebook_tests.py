from pathlib import Path

import pytest
from testbook import testbook


@pytest.mark.parametrize(
    "notebook",
    [
        "HTGR",
        "IAEA3D",
        pytest.param(
            "PWRMOX",
            marks=pytest.mark.skip("Currently unsupported due to missing files"),
        ),
        "PWR_rod_bundle",
    ],
)
def test_notebook(notebook):
    d = Path("../docs/source")
    with testbook(str(d / f"{notebook}.ipynb"), execute=True) as tb:
        assert tb.ref("k")
