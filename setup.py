from setuptools import setup, find_packages
import glob

with open('envs/requirements.txt') as file_open:
     requirements = file_open.read().splitlines()

setup(
    name="gplas",
    setup_requires=[
        "setuptools>=38.6.0",
        "setuptools_scm",
        "setuptools_scm_git_archive",
    ],
    use_scm_version={"write_to":"gplas/version.py"},
    #version="1.1.2-beta",
    scripts=[script for script in glob.glob("gplas/scripts/m_*.py")],
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    package_data={'': ['gplas/*']},
    entry_points={
        'console_scripts': [
            'gplas = gplas.gplas:main',
            #'start = gplas.__main__:start',
            #'dostart = gplas:start'
            ],
    }

)


