# Contributing

````{note} TLDR:
Create a conda environment with:
```{code-block} shell
conda env create -f workflows/envs/devel.yml
conda activate GTExSnake
pre-commit install --install-hooks
```
then commit with `cz commit`,
and everything should take care of itself!
````

Welcome, friend!
Open-source software isn't open open-source without the community.
We appreciate your interest and welcome all contributions.
To help keep everything moving smoothly,
we have a few guidelines.

## Bugs

If you think you've found a bug,
let us know [here](HOLDING).
We'll do our best to deal with it ASAP,
but please be patient as we also work many other projects!

## Developing

If you think you can fix one of the bugs,
or would like to submit a new feature,
then let's get coding!

Once you've cloned the repository,
fork it,
and get your development environment set up.
We use [conda](https://docs.conda.io/en/latest/),
[nox](https://nox.thea.codes/en/stable/),
and [pre-commit](https://pre-commit.com/)
to handle environments, testing, and linting.
Between them,
they make sure that all checks run in isolated environments.
Please make sure you activate them before making any changes,
using the commands below:

```{code-block} shell
conda env create -f workflows/envs/devel.yml
conda activate GTExSnake
pre-commit install --install-hooks
```

Now,
you don't even have to think about linting or testing.
When you commit a change,
pre-commit will automatically run [black](https://github.com/psf/black),
[isort](https://pycqa.github.io/isort/),
[mypy](https://mypy.readthedocs.io/en/stable/index.html),
and a suite of [flake8](https://flake8.pycqa.org/en/latest/)-based linters.
When you push a change,
github actions will trigger a more robust suit using Nox,
including security checks,
testing with [pytest](https://docs.pytest.org/en/6.2.x/),
and automatic doc building with [sphinx](https://www.sphinx-doc.org/en/master/).

### Commits

We use a GitHub action for
[python-semantic-release](https://github.com/relekang/python-semantic-release)
to manage our version numbers.
This automatically parses your commit messages to determine if a new release is nesessary,
and, if so, what kind (ie. major, minor, patch).
So your commit message is very important!

But we also don't want you stressing about how to format your commit message.
To that end,
we use a python implementation of
[commitizen](https://commitizen-tools.github.io/commitizen/index.html)
to handle that for you!
Just commit using `cz commit` instead of `git commit`,
and enjoy the magic!

### Documentation and Testing

Speaking of documentation and testing -
if you add new code,
please add documentation and tests for it as well.
We use [napoleon numpy](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html)
for our docstrings.
Include the docstring with the code,
and sphinx will handle everything else!

## Review

Once your happy with your code,
open a pull-request,
and we will reveiw ASAP.
If you pull-request is not passing on github actions,
or something else confuses us
(we are, after all, only human!),
we might ask for some small changes.
Once everything is looking good,
we will merges in your changes!

## From the Command-Line

```{note} Optional fun!
None of this section is required, but we find it useful and hope you do, too!
```

If you haven't heard of it already,
give a peek to to [gh](https://github.com/cli/cli),
GitHub's official CLI.
It allows to manage all of the above steps from the command-line,
from forking,
to raising issues,
and checking on the status of your pull request.
Not a necessity,
but for you terminal warriors out there,
it just might help!
