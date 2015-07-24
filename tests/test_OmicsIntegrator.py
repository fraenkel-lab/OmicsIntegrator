#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_OmicsIntegrator
----------------------------------

Tests for `OmicsIntegrator`: what tests should we run to ensure that future updates
to garnet and forest do not break its usability.  Current use cases:

1-garnet.py: no expression data
2-garnet.py: with expression data
3-forest.py: with garnet data and protein data
4-forest.py: with only protein data
"""

import unittest

from OmicsIntegrator import garnet-forest


class TestOmicsIntegrator(unittest.TestCase):

    def setUp(self):
        pass

    def test_something(self):
        pass

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
