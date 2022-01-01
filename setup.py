"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/tutorials/packaging-projects/
"""
import re
import os
import setuptools
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

# Get version number from __init__.py
here = os.path.abspath(os.path.dirname(__file__))
regex = "(?<=__version__..\s)\S+"
with open(os.path.join(here,'ISRpy/__init__.py'),'r', encoding='utf-8') as f:
    text = f.read()
match = re.findall(regex,text)
version = match[0].strip("'")


setuptools.setup(
    name="ISRpy",
    version=version,
    author="Erhan Kudeki",
    author_email="erhan@illinois.edu",
    description="Radar Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/radars-eceillinois/isrpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy","scipy"],
    include_package_data=True,
    package_data={"ISRpy": ["igrfdata/*.txt",
                       "igrfdata/*.DAT"],
    },
    python_requires=">=3.4",
)
