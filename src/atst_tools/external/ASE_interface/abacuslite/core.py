"""Forward the historical calculator module to the installed implementation."""

from importlib import import_module
import sys

sys.modules[__name__] = import_module("abacuslite.core")
