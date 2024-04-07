#!/usr/bin/env python
from MODApy.cfg import (
    configuration,
)  # Fix for sphinx - have to research why this is needed
from MODApy.version import __version__


"""
A Package to parse and analyze multiple omics data
"""
name = "MODApy"

__version__ = __version__
configuration = configuration
