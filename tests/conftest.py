# content of conftest.py
# See examples in https://pytest.org/latest/example/simple.html and
# https://pytest.org/latest/fixture.html#fixture-function
import pytest

def pytest_addoption(parser):
    # Do not set a default value
    parser.addoption('--msgpath',dest='msgsteiner',type='string',\
        help='Full path to the msgsteiner dependency, including the executable name.')

@pytest.fixture
def msgsteiner(request):
    return request.config.getoption("--msgpath")
