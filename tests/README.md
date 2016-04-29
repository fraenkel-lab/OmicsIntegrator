===============================
Omics Integrator Testing
===============================

This directory contains Forest unit tests.  The primary test`test_integration.py`
runs Forest on example data and compares the outputs to a stored set of standard outputs.
If the example data files have not been downloaded, the integration test will fail
but the other tests can still be run.

To run the tests use the command
```
py.test --msgpath=<msgsteiner_path>
```
where **msgsteiner_path** is the path to the `msgsteiner` executable, including the executable name.

It is also possible to run these tests via `setup.py` from the `OmicsIntegrator` directory
```
python setup.py test -a "--msgpath=<msgsteiner_path>"
```
which passes the `msgpath` argument to the tests.

These tests require the [pytest package](https://pytest.org/latest/getting-started.html).