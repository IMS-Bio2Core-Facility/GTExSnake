# -*- coding: utf-8 -*-
"""Nox session configuration."""
import nox
from nox.sessions import Session

PACKAGE: str = "workflow/scripts"
LOCATIONS: list[str] = [
    PACKAGE,
    "noxfile.py",
]
VERSIONS: list[str] = ["3.9"]
PIP_PARAMS: list[str] = ["-r"]
CONDA_PARAMS: list[str] = ["-c", "bioconda", "-c", "conda-forge", "--file"]

nox.options.stop_on_first_error = True
nox.options.default_venv_backend = "virtualenv"
nox.options.reuse_existing_virtualenvs = False
nox.options.sessions = ["form", "lint", "type", "security", "doc_tests"]


@nox.session(python=VERSIONS[0])
def form(session: Session) -> None:
    """Format code with isort and black."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "environments/form.txt")
    session.run("isort", *args)
    session.run("black", *args)


@nox.session(python=VERSIONS)
def lint(session: Session) -> None:
    """Lint files with flake8."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "environments/lint.txt")
    session.run("pflake8", *args)


@nox.session(python=VERSIONS)
def type(session: Session) -> None:
    """Type check files with mypy."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "environments/type.txt")
    session.run(
        "mypy",
        "--ignore-missing-imports",
        "--disable-error-code",
        "name-defined",
        *args
    )


# The snakemake install is complex, and runs better in conda
@nox.session(python=VERSIONS, venv_backend="conda")
def security(session: Session) -> None:
    """Check security safety."""
    args = session.posargs or []
    session.conda_install(*CONDA_PARAMS, "environments/security.txt")
    session.run("safety", "check", "--bare", *args)


@nox.session(python=VERSIONS)
def doc_tests(session: Session) -> None:
    """Test the docs.

    As xdoctest does not seem to detect if multiple sources are passed,
    session run must be called over each source manually.

    Parameters
    ----------
    session : Session
        nox session

    """
    args = session.posargs or []
    command = [
        "python",
        "-m",
        "xdoctest",
        "--verbose",
        "2",
        "--report",
        "cdiff",
        "--nocolor",
    ]
    session.install(*PIP_PARAMS, "environments/doc_tests.txt")
    for x in LOCATIONS:
        session.run(*command, x, *args)


@nox.session(python="3.9")
def doc_build(session: Session) -> None:
    """Build the docs."""
    session.install(*PIP_PARAMS, "environments/doc_build.txt")
    session.run("sphinx-build", "docs", "docs/_build")
