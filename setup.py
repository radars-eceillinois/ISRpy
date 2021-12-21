"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/tutorials/packaging-projects/
"""

import setuptools
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ISRpy",
    version="0.0.12",
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
