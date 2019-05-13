# ---------------------------------------------------------------------
#
# Copyright (C) 2016 - 2018 by the deal.II authors
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
import math
from PyDealII.Debug import *


class TestPointWrapper(unittest.TestCase):

    def test_2d_point(self):
        p1 = Point([0., 1.])
        self.assertEqual(p1.x, 0.)
        self.assertEqual(p1.y, 1.)
        p1.x = 1.
        p1.y = 2.
        self.assertEqual(p1.x, 1.)
        self.assertEqual(p1.y, 2.)
        p2 = Point([0., 2.])
        self.assertEqual(p1.distance(p2), 1.)
        self.assertEqual(p2.norm(), 2.)
        self.assertEqual(p2.norm_square(), 4.)
        self.assertEqual(p1 != p2, True)
        self.assertEqual(p1 == p2, False)
        self.assertEqual(p1*p2, 4.)
        p3 = p1 + p2
        self.assertEqual(p3.x, p1.x + p2.x)
        self.assertEqual(p3.y, p1.y + p2.y)
        p3 = p1 - p2
        self.assertEqual(p3.x, p1.x - p2.x)
        self.assertEqual(p3.y, p1.y - p2.y)
        p3 = -p2
        self.assertEqual(p3.x, -p2.x)
        self.assertEqual(p3.y, -p2.y)
        p3 = p2 / 2.
        self.assertEqual(p3.x, p2.x / 2.)
        self.assertEqual(p3.y, p2.y / 2.)
        p3 = p2 * 2.
        self.assertEqual(p3.x, p2.x * 2.)
        self.assertEqual(p3.y, p2.y * 2.)
        p2 += p1
        self.assertEqual(p2.x, 1.)
        self.assertEqual(p2.y, 4.)
        p2 -= p1
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 2.)
        p2 /= 2.
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 1.)
        p2 *= 2.
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 2.)

    def test_3d_point(self):
        p1 = Point([0., 1., 2.])
        self.assertEqual(p1.x, 0.)
        self.assertEqual(p1.y, 1.)
        self.assertEqual(p1.z, 2.)
        p1.x = 1.
        p1.y = 2.
        p1.z = 3.
        self.assertEqual(p1.x, 1.)
        self.assertEqual(p1.y, 2.)
        self.assertEqual(p1.z, 3.)
        p2 = Point([0., 1., 2.])
        self.assertAlmostEqual(p1.distance(p2), math.sqrt(3))
        self.assertAlmostEqual(p2.norm(), math.sqrt(5))
        self.assertEqual(p2.norm_square(), 5.)
        self.assertEqual(p1 != p2, True)
        self.assertEqual(p1 == p2, False)
        self.assertEqual(p1*p2, 8)
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
        p3 = p2 / 2.
        self.assertEqual(p3.x, p2.x / 2.)
        self.assertEqual(p3.y, p2.y / 2.)
        self.assertEqual(p3.z, p2.z / 2.)
        p3 = p2 * 2.
        self.assertEqual(p3.x, p2.x * 2.)
        self.assertEqual(p3.y, p2.y * 2.)
        self.assertEqual(p3.z, p2.z * 2.)
        p2 += p1
        self.assertEqual(p2.x, 1.)
        self.assertEqual(p2.y, 3.)
        self.assertEqual(p2.z, 5.)
        p2 -= p1
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 1.)
        self.assertEqual(p2.z, 2.)
        p2 /= 2.
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 0.5)
        self.assertEqual(p2.z, 1.)
        p2 *= 2.
        self.assertEqual(p2.x, 0.)
        self.assertEqual(p2.y, 1.)
        self.assertEqual(p2.z, 2.)


if __name__ == '__main__':
    unittest.main()
