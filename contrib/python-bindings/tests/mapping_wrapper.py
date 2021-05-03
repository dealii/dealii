# ---------------------------------------------------------------------
#
# Copyright (C) 2016 - 2019 by the deal.II authors
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

class TestMappingWrapperCube(unittest.TestCase):

    def setUp(self):
        self.triangulation = Triangulation('2D')
        self.triangulation.generate_hyper_cube()
        self.triangulation.refine_global(1)

    def test_mapping(self):
        mapping = MappingQGeneric(dim = 2, spacedim = 2, degree = 1)
        p_unit = Point([0.5, 0.5])
        
        for cell in self.triangulation.active_cells():
            p_real = mapping.transform_unit_to_real_cell(cell, p_unit)
            p_unit_back = mapping.transform_real_to_unit_cell(cell, p_real)
            
            self.assertTrue(p_unit.distance(p_unit_back) < 1e-10)

class TestMappingWrapperSphere(unittest.TestCase):

    def setUp(self):
        self.triangulation = Triangulation('2D', '3D')
        p_center = Point([0., 0., 0.])
        self.triangulation.generate_hyper_sphere(p_center)

    def test_mapping(self):
        mapping = MappingQGeneric(dim = 2, spacedim = 3, degree = 4)
        p_unit = Point([0.5, 0.5])
        
        for cell in self.triangulation.active_cells():
            p_real = mapping.transform_unit_to_real_cell(cell, p_unit)
            p_unit_back = mapping.transform_real_to_unit_cell(cell, p_real)
            
            self.assertTrue(p_unit.distance(p_unit_back) < 1e-10)

if __name__ == '__main__':
    unittest.main()
