## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2019 - 2025 by the deal.II authors
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


class TestMappingWrapperCube(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("2D")
        self.triangulation.generate_hyper_cube()
        self.triangulation.refine_global(1)

    def test_mapping(self):
        mapping = MappingQ(dim=2, spacedim=2, degree=1)
        p_unit = Point([0.5, 0.5])

        for cell in self.triangulation.active_cells():
            p_real = mapping.transform_unit_to_real_cell(cell, p_unit)
            p_unit_back = mapping.transform_real_to_unit_cell(cell, p_real)

            self.assertTrue(p_unit.distance(p_unit_back) < 1e-10)


class TestMappingWrapperSphere(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("2D", "3D")
        p_center = Point([0.0, 0.0, 0.0])
        self.triangulation.generate_hyper_sphere(p_center)

    def test_mapping(self):
        mapping = MappingQ(dim=2, spacedim=3, degree=4)
        p_unit = Point([0.5, 0.5])

        for cell in self.triangulation.active_cells():
            p_real = mapping.transform_unit_to_real_cell(cell, p_unit)
            p_unit_back = mapping.transform_real_to_unit_cell(cell, p_real)

            self.assertTrue(p_unit.distance(p_unit_back) < 1e-10)


if __name__ == "__main__":
    unittest.main()
