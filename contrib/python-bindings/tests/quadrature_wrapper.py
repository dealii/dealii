## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2020 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

import unittest

try:
    from PyDealII.Debug import *
except ImportError:
    from PyDealII.Release import *


class TestQuadratureWrapper(unittest.TestCase):
    def test_gauss(self):
        quadrature = Quadrature(dim=3)
        quadrature.create_gauss(n=1)

        w_sum = 0.0
        for weight in quadrature.weights():
            w_sum += weight

        self.assertTrue(abs(1.0 - w_sum) < 1e-10)

    def test_gauss_lobatto(self):
        quadrature = Quadrature(dim=3)
        quadrature.create_gauss_lobatto(n=2)

        w_sum = 0.0
        for weight in quadrature.weights():
            w_sum += weight

        self.assertTrue(abs(1.0 - w_sum) < 1e-10)


if __name__ == "__main__":
    unittest.main()
