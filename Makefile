# Makefile for python packages
#
# Author: Juan Carlos Vázquez

PROJECT_NAME := $(shell python setup.py --name)
PROJECT_VERSION := $(shell python setup.py --version)

default:
	@echo 'Makefile for python packages'
	@echo
	@echo 'Usage:'
	@echo 'make clean	clean all temporary files'
	@echo 'make publish	publish changes to github/PyPI'
	@echo

clean:
	rm -Rf $(PROJECT_NAME).egg-info build dist

publish:
	git add -u
	git commit
	git push
	python3 setup.py sdist bdist_wheel
	twine upload dist/*

publish_test:
	python3 setup.py sdist bdist_wheel
	twine upload --repository testpypi dist/*
