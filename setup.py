import setuptools

from conda_setup import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="dorban",
    version="0.2.0",
    author="Guy Shtotland",
    description="Code for reactor analysis using the diffusion approximation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(exclude=("tests*",)),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
    requirements_yml="conda_requirements.yml"
)
