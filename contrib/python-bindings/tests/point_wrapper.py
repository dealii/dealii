# ---------------------------------------------------------------------
#
# Copyright (C) 2016 by the deal.II authors
#
# This file is part of the deal.II library.
#
# The deal.II library is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the deal.II distribution.
#
# ---------------------------------------------------------------------

import unittest
from PyDealII.Debug import *


class TestPointWrapper(unittest.TestCase):

    def test_2d_point(self):
        point = Point([0., 1.])
        self.assertEqual(point.x, 0.)
        self.assertEqual(point.y, 1.)
        point.x = 1.
        point.y = 2.
        self.assertEqual(point.x, 1.)
        self.assertEqual(point.y, 2.)

    def test_3d_point(self):
        point = Point([0., 1., 2.])
        self.assertEqual(point.x, 0.)
        self.assertEqual(point.y, 1.)
        self.assertEqual(point.z, 2.)
        point.x = 1.
        point.y = 2.
        point.z = 3.
        self.assertEqual(point.x, 1.)
        self.assertEqual(point.y, 2.)
        self.assertEqual(point.z, 3.)


if __name__ == '__main__':
    unittest.main()
