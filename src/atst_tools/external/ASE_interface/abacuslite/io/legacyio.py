"""Forward historical legacy-output imports to installed abacuslite."""

from importlib import import_module
import sys

sys.modules[__name__] = import_module("abacuslite.io.legacyio")
