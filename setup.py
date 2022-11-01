import setuptools
from setuptools.command.install import install


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        print("Configuring MODApy Environment")
        print("Step 1: Verifying installed components")
        print("Step 2: Verifying Annotation Databases")
        print("Step 3: Verifying Reference files")
        print('Done! To run MODApy, simply run "MODApy launcher" from the console')
        # PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION
        install.run(self)


with open("README.md", "r") as rd:
    long_description = rd.read()

version = {}
with open("./MODApy/version.py") as v:
    exec(v.read(), version)

setuptools.setup(
    name="MODApy",
    version=version["__version__"],
    author="Juan Carlos Vázquez, Elmer A. Fernández",
    author_email="juancgvazquez@gmail.com",
    description="Package to perform several analysis on Multi-Omics Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juancgvazquez/MODApy",
    packages=setuptools.find_packages(),
    package_data={"MODApy": ["www/*"]},
    install_requires=[
        "pandas>=1.3.0",
        "numpy",
        "configparser",
        "argparse",
        "Cython",
        "cyvcf2==0.9.0",
        "xlrd",
        "openpyxl",
        "XlsxWriter",
        "matplotlib",
        "matplotlib-venn",
        "xmltodict",
        "pyyaml",
        "requests",
        "tqdm",
        "fastapi",
        "rq",
        "uvicorn",
    ],
    entry_points={
        # Command line scripts
        "console_scripts": ["MODApy=MODApy.cmd_line:main"]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
