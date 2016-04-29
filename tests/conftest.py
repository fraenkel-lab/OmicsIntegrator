# content of conftest.py
import pytest

def pytest_addoption(parser):
    parser.addoption('--msgpath',dest='msgsteiner',type='string',help='Full path to the msgsteiner dependency, including the executable name.')

@pytest.fixture
def msgsteiner(request):
    return request.config.getoption("--msgpath")
