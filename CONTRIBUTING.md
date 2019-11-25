# Contributing to cellmlmanip

## Testing

To run tests, just run
```sh
pytest
```

### Automated tests

Automated testing is run on linux using Travis CI.

### Writing tests

Tests should be grouped and named in a manner that makes it easy to find the test for an appropriate piece of code.
As in the main code, docstrings and comments should be added to clarify the intention of a piece of code (allowing readers to check if it does what it's supposed to do).

## Coding style

We aim to be [PEP8](https://www.python.org/dev/peps/pep-0008/) compliant, with a slightly more relaxed attitude towards line length. Lines up to 100 chars are allowed if it means better readability.

### flake8

`flake8` carries out PEP8 and other coding style checks. 
It can be run by typing in `flake8` and the output of this command should be kept clean.

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


## Documentation

Documentation is generated from the docstrings in the code, and files in the ``docs`` directory.
Both are written using [reStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html).
This is then compiled into HTML (or other formats) using [Sphinx](http://www.sphinx-doc.org/en/stable/).
We use [readthedocs.org](https://readthedocs.org) to do this automatically, with the result visible at [https://cellmlmanip.readthedocs.io](https://cellmlmanip.readthedocs.io).

### Building the documentation locally

Install with `pip install -e .[docs]`,  navigate to the `docs` directory and build with `make html`.

### Docstring guidelines

Ideally, docstrings start with a single line comment explaining what the method does.
This is followed by a blank line, and an optional longer description of how the method works.
For complex methods, arguments and return types can be clarified using the ``:param number: <description>`` syntax.
An example is given below:

```
def make_tea(self, n_cups, add_milk=False):
    """
    Brews and returns ``n_cups`` of tea.
    
    Tea is brewed using Barry's algorithm using a specially selected blend.
    See also :meth:`cellmlmanip.eat_biscuit()`.
    
    :param n_cups: The (integer) number of cups of tea to brew.
    :param add_milk: Set to ``True`` to add milk to all returned cups of tea.
    :return: An iterator over :class:`CupOfTea` objects.
    """
```
