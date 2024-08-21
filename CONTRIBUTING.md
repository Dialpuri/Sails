# Contributing

Contributions are welcome, please open a pull request with any changes and add Dialpuri as a reviewer.

## Development

To build Sails upon every Python import, which will be useful in development, run:

    pip install --no-build-isolation --config-settings=editable.rebuild=true -Cbuild-dir=build -ve .

With CLion, to get proper intellisense, load the CMake Project with this additional setting

    -Dnanobind_DIR=.venv/lib/python3.10/site-packages/nanobind/cmake


[//]: # (## Testing)

[//]: # (Any changes must pass the tests defined in `package/tests`. Test can be ran using `pytest` with: )

[//]: # ()
[//]: # (    pytest package/tests --unit --runslow -v)

[//]: # ()
