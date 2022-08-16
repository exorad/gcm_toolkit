"""
General GCMT writer
"""
import os

from gcm_toolkit import GCMT
from gcm_toolkit.core.writer import Writer, write_message


def test_writer_with_file():
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    outputfile = "testfile.txt"
    GCMT(write=outputfile)
    assert os.path.exists(outputfile)
    os.remove(outputfile)


def test_writer_wrong():
    """Writer should be disabled in that case"""
    GCMT(write=2)
    assert not Writer.on


def test_writer_linebreak():
    """Writer should do automatic linebreaks"""
    GCMT(write="on")
    write_message(100 * "test")
    write_message(100 * " ")
