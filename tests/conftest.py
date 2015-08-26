# content of conftest.py
import pytest

def pytest_addoption(parser):
    parser.addoption('--msgsteiner',dest='msgsteiner',type='string',help='Path to msgsteiner9 code, be sure to include!')

@pytest.fixture
def msgsteiner(request):
    return request.config.getoption("--msgsteiner")
