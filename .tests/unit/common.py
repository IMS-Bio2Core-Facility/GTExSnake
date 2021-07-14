# -*- coding: utf-8 -*-
"""
Common code for unit testing of rules generated with Snakemake 6.4.0.
"""

import hashlib
import os
import subprocess as sp
import sys
from difflib import unified_diff
from pathlib import Path


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if (
                    str(f).startswith(".snakemake")
                    or str(f).endswith(".txt")
                    or str(f).endswith(".log")
                    or str(f).startswith("config")
                    or str(f).startswith(
                        "resources"
                    )  # these turn up in integration tests
                ):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])


def file_diff(self, gen_file, ex_file):
    # with block for file opening prevents easy comprehension
    contents = []
    for file in [gen_file, ex_file]:
        with open(file, "r") as f:
            contents.append(f.readlines())
    sys.stderr.writelines(
        unified_diff(*contents, fromfile="Generated File", tofile="Expected File")
    )
