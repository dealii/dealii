## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2016 - 2025 by the deal.II authors
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
import math

try:
    from PyDealII.Debug import *
except ImportError:
    from PyDealII.Release import *


class TestPointWrapper(unittest.TestCase):
    def test_2d_point(self):
        p1 = Point([0.0, 1.0])
        self.assertEqual(p1.x, 0.0)
        self.assertEqual(p1.y, 1.0)
        p1.x = 1.0
        p1.y = 2.0
        self.assertEqual(p1.x, 1.0)
        self.assertEqual(p1.y, 2.0)
        p2 = Point([0.0, 2.0])
        self.assertEqual(p1.distance(p2), 1.0)
        self.assertEqual(p2.norm(), 2.0)
        self.assertEqual(p2.norm_square(), 4.0)
        self.assertEqual(p1 != p2, True)
        self.assertEqual(p1 == p2, False)
        self.assertEqual(p1 * p2, 4.0)
        p3 = p1 + p2
        self.assertEqual(p3.x, p1.x + p2.x)
        self.assertEqual(p3.y, p1.y + p2.y)
        p3 = p1 - p2
        self.assertEqual(p3.x, p1.x - p2.x)
        self.assertEqual(p3.y, p1.y - p2.y)
        p3 = -p2
        self.assertEqual(p3.x, -p2.x)
        self.assertEqual(p3.y, -p2.y)
        p3 = p2 / 2.0
        self.assertEqual(p3.x, p2.x / 2.0)
        self.assertEqual(p3.y, p2.y / 2.0)
        p3 = p2 * 2.0
        self.assertEqual(p3.x, p2.x * 2.0)
        self.assertEqual(p3.y, p2.y * 2.0)
        p2 += p1
        self.assertEqual(p2.x, 1.0)
        self.assertEqual(p2.y, 4.0)
        p2 -= p1
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 2.0)
        p2 /= 2.0
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 1.0)
        p2 *= 2.0
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 2.0)

    def test_3d_point(self):
        p1 = Point([0.0, 1.0, 2.0])
        self.assertEqual(p1.x, 0.0)
        self.assertEqual(p1.y, 1.0)
        self.assertEqual(p1.z, 2.0)
        p1.x = 1.0
        p1.y = 2.0
        p1.z = 3.0
        self.assertEqual(p1.x, 1.0)
        self.assertEqual(p1.y, 2.0)
        self.assertEqual(p1.z, 3.0)
        p2 = Point([0.0, 1.0, 2.0])
        self.assertAlmostEqual(p1.distance(p2), math.sqrt(3))
        self.assertAlmostEqual(p2.norm(), math.sqrt(5))
        self.assertEqual(p2.norm_square(), 5.0)
        self.assertEqual(p1 != p2, True)
        self.assertEqual(p1 == p2, False)
        self.assertEqual(p1 * p2, 8)
        dim = 3
        p3 = p1 + p2
        self.assertEqual(p3.x, p1.x + p2.x)
        self.assertEqual(p3.y, p1.y + p2.y)
        self.assertEqual(p3.z, p1.z + p2.z)
        p3 = p1 - p2
        self.assertEqual(p3.x, p1.x - p2.x)
        self.assertEqual(p3.y, p1.y - p2.y)
        self.assertEqual(p3.z, p1.z - p2.z)
        p3 = -p2
        self.assertEqual(p3.x, -p2.x)
        self.assertEqual(p3.y, -p2.y)
        self.assertEqual(p3.z, -p2.z)
        p3 = p2 / 2.0
        self.assertEqual(p3.x, p2.x / 2.0)
        self.assertEqual(p3.y, p2.y / 2.0)
        self.assertEqual(p3.z, p2.z / 2.0)
        p3 = p2 * 2.0
        self.assertEqual(p3.x, p2.x * 2.0)
        self.assertEqual(p3.y, p2.y * 2.0)
        self.assertEqual(p3.z, p2.z * 2.0)
        p2 += p1
        self.assertEqual(p2.x, 1.0)
        self.assertEqual(p2.y, 3.0)
        self.assertEqual(p2.z, 5.0)
        p2 -= p1
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 1.0)
        self.assertEqual(p2.z, 2.0)
        p2 /= 2.0
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 0.5)
        self.assertEqual(p2.z, 1.0)
        p2 *= 2.0
        self.assertEqual(p2.x, 0.0)
        self.assertEqual(p2.y, 1.0)
        self.assertEqual(p2.z, 2.0)


if __name__ == "__main__":
    unittest.main()
