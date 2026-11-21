#!/usr/bin/env python3
"""End-to-end regression test for mo_irreps."""

import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "tools" / "misc" / "mo_irreps.py"
DATA = Path(__file__).with_name("methane")
SPEC = importlib.util.spec_from_file_location("mo_irreps", SCRIPT)
MO_IRREPS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MO_IRREPS
SPEC.loader.exec_module(MO_IRREPS)


class MoIrrepsRegressionTest(unittest.TestCase):
    """Run one complete molecular symmetry analysis."""

    def test_tetrahedral_methane(self):
        """Recover occupied and virtual A1 and T2 methane manifolds."""
        output = io.StringIO()
        arguments = [
            "--eigenvec",
            str(DATA / "eigenvec.out"),
            "--geometry",
            str(DATA / "ref.xyz"),
            "--band",
            str(DATA / "band.out"),
            "--overlap",
            str(DATA / "oversqr.dat"),
            "--exc-input",
            str(DATA / "missing-exc.dat"),
        ]
        with contextlib.redirect_stdout(output):
            self.assertEqual(MO_IRREPS.main(arguments), 0)

        lines = [
            line.split()
            for line in output.getvalue().splitlines()
            if line and not line.startswith("#")
        ]
        result = [
            (fields[2], float(fields[3]), float(fields[4]), fields[5])
            for fields in lines
        ]
        self.assertEqual(
            result,
            [
                ("1", -15.633, 2.0, "A1"),
                ("2,3,4", -9.14, 6.0, "T2"),
                ("5,6,7", 9.515, 0.0, "T2"),
                ("8", 16.118, 0.0, "A1"),
            ],
        )


if __name__ == "__main__":
    unittest.main()
