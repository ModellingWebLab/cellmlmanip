[![Build Status](https://travis-ci.org/ModellingWebLab/cellmlmanip.svg?branch=develop)](https://travis-ci.org/ModellingWebLab/cellmlmanip)

# cellmlmanip
CellML loading and model equation manipulation

## Developer installation

Set up a virtual environment, and install requirements:
```sh
conda create -n cellmlmanip python=3.6
source activate cellmlmanip
pip install -r requirements/dev.txt
pip install -e .
```

## Testing

To run tests, just run
```sh
pytest
```

## Coding style

We aim to be [PEP8](https://www.python.org/dev/peps/pep-0008/) compliant, with a slightly more relaxed attitude towards line length. Lines up to 100 chars are allowed if it means better readability.

### flake8

`flake8` carries out PEP8 and other coding style checks. It can be run by typing in `flake8` and the output of this command should be kept clean.

It is recommended to install the git pre-commit hook:

```sh
flake8 --install-hook=git
```

By default, this will warn of code style violations before you commit, but will still allow the commit. You can set it to strict mode:

```sh
git config --bool flake8.strict true
```

This will stop the commit from happening if flake8 checks do not pass.


### isort

`isort` keeps imports in the [PEP8 recommended order](https://www.python.org/dev/peps/pep-0008/#id23) (flake8 / pep8 tools don't handle this). It can be run manually by typing in `isort -rc .` in `cellmlmanip/` and this should generate no errors.

There is a git pre-commit hook for this too, the process is a little more manual than for flake8:

https://github.com/timothycrosley/isort#git-hook
