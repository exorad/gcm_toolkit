"""
Common functions for gcm_toolkit testing
"""

import fnmatch
import hashlib
import os
import tarfile
import tempfile
import urllib.request as req
from contextlib import contextmanager

import py
import pytest
from xmitgcm.file_utils import clear_cache


@contextmanager
def hide_file(origdir, *basenames):
    """Temporarily hide files within the context."""
    # make everything a py.path.local
    tmpdir = py.path.local(tempfile.mkdtemp())
    origdir = py.path.local(origdir)
    oldpaths = [origdir.join(basename) for basename in basenames]
    newpaths = [tmpdir.join(basename) for basename in basenames]

    # move the files
    for oldpath, newpath in zip(oldpaths, newpaths):
        oldpath.rename(newpath)
    # clear the cache if it exists
    clear_cache()
    try:
        yield str(tmpdir)
    finally:
        # move them back
        for oldpath, newpath in zip(oldpaths, newpaths):
            newpath.rename(oldpath)


DLROOT = "https://ndownloader.figshare.com/files/"

# parameterized fixture are complicated
# http://docs.pytest.org/en/latest/fixture.html#fixture-parametrize

# dictionary of archived experiments and some expected properties
_experiments = {
    "HD2_test": {
        "dlink": DLROOT + "36234516",
        "md5": "20a49edac60f905cffdd1300916e978c",
        "gcm": "MITgcm",
        "rel_data_dir": "{}/run",
        "p_domain": [1e-5, 650],  # expected pressure domain in bar
        "times": [12000],  # expected timestamps in days
        "area_key": "area_c",
        "v_data": "V",
    },
}


def setup_mds_dir(tmpdir_factory, request, dbs):
    """Helper function for setting up test cases."""
    expt_name = request.param
    expected_results = dbs[expt_name]
    target_dir = str(tmpdir_factory.mktemp("mdsdata"))
    try:
        # user-defined directory for test datasets
        data_dir = os.environ["GCMTOOLS_TESTDATA"]
    except KeyError:
        # default to HOME/.xmitgcm-test-data/
        data_dir = os.environ["HOME"] + "/.gcmtools-test-data"
    # create the directory if it doesn't exixt
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    datafile = os.path.join(data_dir, expt_name + ".tar.gz")
    # download if does not exist locally
    if not os.path.exists(datafile):
        print("File does not exist locally, downloading...")
        download_archive(expected_results["dlink"], datafile)
        localmd5 = file_md5_checksum(datafile)
        if localmd5 != expected_results["md5"]:
            os.remove(datafile)
            msg = """
            MD5 checksum does not match, try downloading dataset again.
            """
            raise IOError(msg)

    return untar(data_dir, expt_name, target_dir), expected_results


def download_archive(url, filename):
    """download file from url into datafile

    PARAMETERS:

    url: str
        url to retrieve
    filename: str
        file to save on disk
    """

    req.urlretrieve(url, filename)


def untar(data_dir, basename, target_dir):
    """Unzip a tar file into the target directory. Return path to unzipped
    directory."""
    datafile = os.path.join(data_dir, basename + ".tar.gz")
    if not os.path.exists(datafile):
        raise IOError("Could not find data file " + datafile)
    with tarfile.open(datafile) as tar:
        tar.extractall(target_dir)
    # subdirectory where file should have been untarred.
    # assumes the directory is the same name as the tar file itself.
    # e.g. testdata.tar.gz --> testdata/
    fulldir = os.path.join(target_dir, basename)
    if not os.path.exists(fulldir):
        raise IOError(r"Could not find tar file output dir " + fulldir)
    # the actual data lives in a file called testdata
    # clean up ugly weird hidden files that mac-os sometimes puts in the archive
    # https://unix.stackexchange.com/questions/9665/create-tar-archive-of-a-directory-except-for-hidden-files
    # https://superuser.com/questions/259703/get-mac-tar-to-stop-putting-filenames-in-tar-archives
    bad_files = [f for f in os.listdir(fulldir) if fnmatch.fnmatch(f, "._*")]
    for fil in bad_files:
        os.remove(os.path.join(fulldir, fil))

    return fulldir


def file_md5_checksum(fname):
    """update md5 file"""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as fil:
        hash_md5.update(fil.read())
    return hash_md5.hexdigest()


# find the tar archive in the test directory
# http://stackoverflow.com/questions/29627341/pytest-where-to-store-expected-data
@pytest.fixture(scope="module", params=_experiments.keys())
def all_raw_testdata(tmpdir_factory, request):
    """setup test data"""
    return setup_mds_dir(tmpdir_factory, request, _experiments)


@pytest.fixture(scope="module", params=["HD2_test"])
def exorad_testdata(tmpdir_factory, request):
    """set up test data"""
    return setup_mds_dir(tmpdir_factory, request, _experiments)
