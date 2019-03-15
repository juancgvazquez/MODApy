import setuptools

with open("README.md", "r") as rd:
    long_description = rd.read()

version = {}
with open('./MODApy/version.py') as v:
    exec(v.read(), version)
setuptools.setup(
    name="MODApy",
    version=version['__version__'],
    author='Juan Carlos Vázquez, Elmer A. Fernández',
    author_email='juancgvazquez@gmail.com',
    description='Package to perform several analysis on Multi-Omics Data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juancgvazquez/MODApy",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'configparser',
        'argparse',
        'Cython',
        'cyvcf2==0.9.0',
        'xlrd',
        'openpyxl',
        'XlsxWriter',
        'matplotlib',
        'matplotlib-venn',
        'xmltodict',
        'pyyaml',
        'requests',
        'tqdm'
    ],
    entry_points={
        # Command line scripts
        'console_scripts': ['MODApy=MODApy.cmd_line:main']
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
