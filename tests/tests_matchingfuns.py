#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the matchingfuns module.
"""

import pytest
import matchingfuns as mat

def test_matchingfuns_errors():
    """
    Testing improper inputs to matchingfuns functions.
    """
    
    with pytest.raises(AssertionError):
        mat.min_switch_ind_finder('foo',[0,2],[0,2],[0,2])