from setuptools import setup, find_packages

with open("requirements.txt") as file_open:
     requirements = file_open.read().splitlines()

with open("README.md") as file_open:
     README = file_open.read()

setup(
    name="gplas",
    setup_requires=[
        "setuptools>=38.6.0",
        "setuptools_scm",
        "setuptools_scm_git_archive",
    ],
    use_scm_version={"version_file":"gplas/version.py"},
    #version="1.0.0",
    description="Binning plasmid-predicted contigs using short-read graphs",
    long_description=README,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    entry_points={
        'console_scripts': ["gplas = gplas.gplas:main"],
    }

)
