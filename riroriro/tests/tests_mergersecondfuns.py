#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the mergersecondfuns module.
"""

import pytest
import riroriro.mergersecondfuns as me2

def test_mergersecondfuns_errors():
    """
    Testing improper inputs to mergersecondfuns functions.
    """
    
    with pytest.raises(AssertionError):
        me2.merger_phase_calculation(1.4,1,[0,2],[0.2])
        
    print('test_mergersecondfuns_errors completed successfully')