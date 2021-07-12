# -*- coding: utf-8 -*-
import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_biomart():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/biomart/data")
        expected_path = PurePosixPath(".tests/unit/biomart/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "results/biomart/ASCL1_message.csv results/biomart/BSX_message.csv results/biomart/CXXC5_message.csv results/biomart/DLX1_message.csv results/biomart/ESR1_message.csv",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/biomart/ASCL1_message.csv results/biomart/BSX_message.csv results/biomart/CXXC5_message.csv results/biomart/DLX1_message.csv results/biomart/ESR1_message.csv",
                "-F",
                "-j1",
                "--keep-target-files",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
