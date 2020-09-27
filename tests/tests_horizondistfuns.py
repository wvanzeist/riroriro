#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the horizondistfuns module.
"""

import numpy as np
import pytest
import horizondistfuns as hor

def test_horizondistfuns_errors():
    """
    Testing improper inputs to horizondistfuns functions.
    """
    
    with pytest.raises(AssertionError):
        hor.horizon_distance_calculation(np.empty((3)),np.empty((3)),\
                                         np.empty((3)),'quad')