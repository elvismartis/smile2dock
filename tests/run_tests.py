#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test runner for SMILES-to-3D converters
"""

import unittest
import sys
import os

def run_tests():
    """Run all tests and return the number of failures"""
    # Add parent directory to path so we can import the modules
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    # Find and run all tests
    loader = unittest.TestLoader()
    start_dir = os.path.dirname(os.path.abspath(__file__))
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return len(result.failures) + len(result.errors)

if __name__ == '__main__':
    sys.exit(run_tests())