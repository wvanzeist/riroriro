# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = []
"""
from .example_mod import *   # noqa
# Then you can be explicit to control what ends up in the namespace,
__all__ += ['do_primes']   # noqa
# or you can keep everything from the subpackage with the following instead
# __all__ += example_mod.__all__
"""