#
# Tests the versioning of cellmlmanip
#
# import pytest
import logging


# Show more logging output
logging.getLogger().setLevel(logging.INFO)


def test_version():
    # Test the version() method
    import cellmlmanip

    version = cellmlmanip.version()
    assert isinstance(version, tuple)
    assert len(version) == 3
    assert isinstance(version[0], int)
    assert isinstance(version[1], int)
    assert isinstance(version[2], int)
    assert version[0] >= 0
    assert version[1] >= 0
    assert version[2] >= 0

    version = cellmlmanip.version(formatted=True)
    assert isinstance(version, str)
    assert len(version) >= 1
    assert version.startswith('cellmlmanip ')
