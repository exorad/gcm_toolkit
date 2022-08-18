""" set up file """
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

SETUP_REQUIRES = ["pytest-runner"]
TESTS_REQUIRE = ["pytest >= 4.0", "coverage"]
INSTALL_REQUIRES = [
    "scipy>=1.7.0",
    "numpy",
    "f90nml",
    "astropy",
    "xarray",
    "pyyaml",
    "pre-commit",
]

setup(
    name="gcm_toolkit",
    version="v0.1.6",
    packages=find_packages(),
    include_package_data=True,
    scripts=["bin/convert_to_gcmt"],
    url="https://github.com/exorad/gcmt",
    license="MIT",
    author="Aaron David Schneider",
    author_email="aarondavid.schneider@nbi.ku.dk",
    description="postprocessing stuff",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=INSTALL_REQUIRES,
    setup_requires=SETUP_REQUIRES,
    tests_require=TESTS_REQUIRE,
)
