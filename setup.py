import setuptools

with open("README.md", "r") as rd:
	long_description = rd.read()

setuptools.setup(
	name="MoDAPy",
	version='0.0.5dev34',
	author='Juan Carlos Vázquez',
	author_email='juancgvazquez@gmail.com',
	description='Package to perform several analysis on Multi-Omics Data',
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/juancgvazquez/MoDAPy",
	packages=setuptools.find_packages(),
	include_package_data=True,
	install_requires=[
		'pandas',
		'configparser',
		'argparse',
		'Cython',
		'cyvcf2==0.9.0',
		'xlrd',
		'openpyxl',
		'XlsxWriter',
	],
	scripts=['scripts/MoDAPy'],
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
