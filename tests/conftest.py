import pytest
import os

this_file = os.path.dirname(__file__)

pytest.TEST_FOLDER = os.path.join(this_file, "../testsets")

def pytest_addoption(parser):
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="run slow tests"
    )
