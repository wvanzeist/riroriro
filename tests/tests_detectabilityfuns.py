#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the detectabilityfuns module.
"""

import pytest
import detectabilityfuns as det

def test_detectabilityfuns_errors():
    """
    Testing improper inputs to detectabilityfuns functions.
    """
    
    with pytest.raises(AssertionError):
        det.cdf_generator(123.4)