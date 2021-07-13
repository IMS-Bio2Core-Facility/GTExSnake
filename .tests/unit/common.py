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


class ShaChecker(OutputChecker):
    def _sha256sum(self, file):
        with open(file, "r") as f:
            val = hashlib.sha256("".join([l for l in f.readlines()]).encode("utf8"))
        return val.hexdigest()

    def _compare_sha256(self, generated_file, expected_file):
        assert self._sha256sum(generated_file) == self._sha256sum(
            expected_file
        ), "sha256 checksums do not match"

    def _file_diff(self, gen_file, ex_file):
        # with block for file opening prevents easy comprehension
        contents = []
        for file in [gen_file, ex_file]:
            with open(file, "r") as f:
                contents.append(f.readlines())
        sys.stderr.writelines(
            unified_diff(*contents, fromfile="Generated File", tofile="Expected File")
        )

    def compare_files(self, generated_file, expected_file):
        # Only needed for the gtf file
        if "MANE" in str(generated_file) or "gencode" in str(generated_file):
            try:
                self._compare_sha256(generated_file, expected_file)
            except AssertionError:
                self._file_diff(generated_file, expected_file)
        else:
            sp.check_output(["cmp", generated_file, expected_file])
