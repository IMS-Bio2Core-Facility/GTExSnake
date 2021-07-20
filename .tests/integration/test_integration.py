# -*- coding: utf-8 -*-
import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_integration():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/integration/data")
        expected_path = PurePosixPath(".tests/integration/expected")
        config_path = PurePosixPath(".tests/integration/config")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")

        # dbg
        print(
            "results/process/ASCL1_isoforms.csv results/process/BSX_isoforms.csv",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/process/ASCL1_isoforms.csv",
                "results/process/BSX_isoforms.csv",
                "-j1",
                "--keep-target-files",
                "--use-singularity",
                "--use-conda",
                "--conda-frontend",
                "mamba",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
