'''
Setup script for beastify.
'''
import pathlib
from setuptools import setup
from setuptools import find_packages

from beastify import __VERSION__ as version_string

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="beastify",
    version=version_string,
    description="Partition your alignment into distinct codon positions and non-coding positions",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/andersgs/beastify",
    author="Anders Gonçalves da Silva",
    author_email="andersgs@gmail.com",
    classifiers=[
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=["click", "pandas", "numpy", "biopython"],
    entry_points={
        "console_scripts": [
            "beastify=beastify.__main__:beastify",
        ]
    },
)
