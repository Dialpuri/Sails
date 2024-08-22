# Contributing

Contributions are welcome, please open a pull request with any changes and add Dialpuri as a reviewer.

## Development

### Editable Installs
To build Sails upon every Python import, which will be useful in development, run:

    pip install --no-build-isolation --config-settings=editable.rebuild=true -Cbuild-dir=build -ve .

Sails uses the latest experimental `setuptools` implementation for `scikit-build-core`, and as such editable installs
are not fully supported. To get around this, Sails can use the traditional `scikit-build-core` build system
for development without much issue, this is located in the `_pyproject.toml` file in the `package` directory. To swap
between the `setuptools` backend (files `setup.py` and `pyproject.toml`) run the `swap.py` script in the `package`
directory. It is imperative that the `setuptools` is used for any deployments so before committing, ensure that
the `test_setuptools.py` test passes. These tests are also run during pull requests.

### Intellisense
With CLion, to get proper intellisense, load the CMake Project with this additional setting

    -Dnanobind_DIR=.venv/lib/python3.X/site-packages/nanobind/cmake

### Formatting

Sails uses `ruff` formatting for Python files, but all files should follow these general rules:
* No trailing whitespace
* Files end with newline
* Double quotes for string literals

These formatting requirements are required before pull requests will be accepted, they are easy to check with

    ruff check package

or fix with

    ruff format

### Pre-Commit  Hooks
It is easy to forget to run formatting checks so it may be a good idea to utilise pre-commit hooks to do so, you can do that in the following way

    pip install pre-commit
    pre-commit install

Following these commands, when items in the `package` directory change and you commit them, the `ruff` formatter will run
and fix any errors, you can then commit the new changes as you would normally.

## Testing

Any changes must pass the tests defined in `package/tests`. Test can be ran using `pytest` with:

    pytest package/tests

Please ensure all tests pass.
