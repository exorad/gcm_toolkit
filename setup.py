from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='GCMtools',
    version='v0.1',
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/exorad/GCMtools',
    license='MIT',
    author='Aaron David Schneider',
    author_email='aarondavid.schneider@nbi.ku.dk',
    description='postprocessing stuff',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "scipy>=1.7.0",
        "numpy",
        "f90nml",
    ]
)
