# ---------------------------------------------------------------------
#
# Copyright (C) 2016 - 2020 by the deal.II authors
#
# This file is part of the deal.II library.
#
# The deal.II library is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE.md at
# the top level directory of deal.II.
#
# ---------------------------------------------------------------------

import unittest
from PyDealII.Debug import *

class TestQuadratureWrapper(unittest.TestCase):

    def test_gauss(self):
        quadrature = Quadrature(dim = 3)
        quadrature.create_gauss(n = 1);
        
        w_sum = 0.
        for weight in quadrature.weights():
            w_sum += weight
            
        self.assertTrue(abs(1. - w_sum) < 1e-10)

    def test_gauss_lobatto(self):
        quadrature = Quadrature(dim = 3)
        quadrature.create_gauss_lobatto(n = 2);
        
        w_sum = 0.
        for weight in quadrature.weights():
            w_sum += weight
            
        self.assertTrue(abs(1. - w_sum) < 1e-10)

if __name__ == '__main__':
    unittest.main()
