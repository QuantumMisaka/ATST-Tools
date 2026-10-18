"""Forward historical K-point helper imports to installed abacuslite."""

from importlib import import_module
import sys

sys.modules[__name__] = import_module("abacuslite.utils.ksampling")
