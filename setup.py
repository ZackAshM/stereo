from os import path
from setuptools import setup

# the stereo version
__version__ = "0.0.1"

# get the absolute path of this project
here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

# the standard setup info
setup(
    name="stereo",
    version=__version__,
    description="A daytime star tracker for the Payload for Ultra-high Energy Observations (PUEO)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rprechelt/stereo",
    author="Remy L. Prechelt",
    author_email="prechelt@hawaii.edu",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    keywords=["star", "tracker", "camera", "astronomy", "balloon", "orbital"],
    packages=["stereo"],
    python_requires=">=3.6*, <4",
    install_requires=[
        "numpy",
        "astropy",
        "matplotlib",
        "uproot",
        "scipy",
        "attrs",
        "timeout_decorator",
        "tqdm",
        "tetra3 @ git+https://github.com/esa/tetra3"
    ],
    extras_require={
        "test": [
            "pytest",
            "black",
            "mypy",
            "isort",
            "coverage",
            "pytest-cov",
            "flake8",
        ],
    },
    scripts=["scripts/stereo"],
    project_urls={},
    include_package_data=True,
)
