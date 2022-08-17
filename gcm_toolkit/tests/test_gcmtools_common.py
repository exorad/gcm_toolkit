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
_experiments_raw = {
    "HD2_test": {
        "dlink": DLROOT + "36234516",
        "md5": "20a49edac60f905cffdd1300916e978c",
        "gcm": "MITgcm",
        "rel_data_dir": "{}/run",
        "p_domain": [1e-5, 650],  # expected pressure domain in bar
        "times": [12000],  # expected timestamps in days
        "iters": [41472000],
        "MMW": 2.3160063003212366,  # in mass of u
        "Rstar": 83692710000.0,  # in cm
        "Tstar": 6092.0,  # in K
        "semimajoraxis": 710141092212.9,  # in cm
        "prt_max": 0.0042081,  # values for petitradtrans test
        "prt_min": 0.00023373,  # values for petitradtrans test
        "prt_mean": 0.00150488,  # values for petitradtrans test
    },
}

_experiments_nc = {
    "HD2_test_nc": {
        "tag": "HD2_test_nc",  # same as the filename (without extension)
        "dlink": DLROOT + "36701961",
        "md5": "73cdd7b938b25a173eeefb4262c6b6db",
        "p_domain": [1e-5, 650],  # expected pressure domain in bar
        "times": [11000, 12000],  # expected timestamps in days
        "iters": [38016000, 41472000],
        "MMW": 2.3160063003212366,  # in mass of u
        "Rstar": 83692710000.0,  # in cm
        "Tstar": 6092.0,  # in K
        "semimajoraxis": 710141092212.9,  # in cm
        "area_key": "area_c",
        "v_data": "V",
        "prt_max": 0.0042081,  # values for petitradtrans test
        "prt_min": 0.00023373,  # values for petitradtrans test
        "prt_mean": 0.00150488,  # values for petitradtrans test
    },
}

_interface_data = {
    "prt_input_data": {
        "dlink": DLROOT + "36704253",
        "md5": "a7a1ecfa6764c0a4b55029da359f42d6",
        "line_species": ["H2O_Exomol_R_1"],
        "rayleigh_species": [],
        "continuum_opacities": [],
        "wlen_bords_micron": [0.3, 30.0],
    },
}


def setup_experiment_dir(tmpdir_factory, request, dbs, format):
    """Helper function for setting up test cases."""
    expt_name = request.param
    expected_results = dbs[expt_name]
    target_dir = str(tmpdir_factory.mktemp("mdsdata"))
    try:
        # user-defined directory for test datasets
        data_dir = os.environ["GCM_TOOLKIT_TESTDATA"]
    except KeyError:
        # default to HOME/.gcm_toolkit-test-data/
        data_dir = os.environ["HOME"] + "/.gcm_toolkit-test-data"
    # create the directory if it doesn't exixt
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    datafile = os.path.join(data_dir, expt_name + format)
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

    if format == ".tar.gz":
        return untar(data_dir, expt_name, target_dir), expected_results
    else:
        return datafile, expected_results


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
@pytest.fixture(scope="module", params=_experiments_nc.keys())
def all_nc_testdata(tmpdir_factory, request):
    """setup test data"""
    return setup_experiment_dir(
        tmpdir_factory, request, _experiments_nc, ".nc"
    )


@pytest.fixture(scope="module", params=_experiments_raw.keys())
def all_raw_testdata(tmpdir_factory, request):
    """setup test data"""
    return setup_experiment_dir(
        tmpdir_factory, request, _experiments_raw, ".tar.gz"
    )


@pytest.fixture(scope="module", params=["HD2_test"])
def exorad_testdata(tmpdir_factory, request):
    """set up test data"""
    return setup_experiment_dir(
        tmpdir_factory, request, _experiments_raw, ".tar.gz"
    )


@pytest.fixture(scope="module", params=["HD2_test_nc"])
def exorad_testdata_nc(tmpdir_factory, request):
    """set up test data"""
    return setup_experiment_dir(
        tmpdir_factory, request, _experiments_nc, ".nc"
    )


@pytest.fixture(scope="module", params=["prt_input_data"])
def petitradtrans_testdata(tmpdir_factory, request):
    """set up test data"""
    return setup_experiment_dir(
        tmpdir_factory, request, _interface_data, ".tar.gz"
    )
