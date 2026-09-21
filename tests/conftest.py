"""Shared test configuration.

Tests import the installed package (``pip install -e .``); nothing here
manipulates ``sys.path``. Tests that need inputs outside the repository --
production four-vectors, decay templates, FONLL variation grids, third-party
builds -- are marked ``external`` and skipped unless the environment provides
them.
"""
import pytest


def pytest_collection_modifyitems(config, items):
    if config.getoption("-m"):
        return
    skip = pytest.mark.skip(reason="needs inputs outside the repository; run with -m external")
    for item in items:
        if "external" in item.keywords:
            item.add_marker(skip)
