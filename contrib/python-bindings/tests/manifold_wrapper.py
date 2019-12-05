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
# The full text of the license can be found in the file LICENSE.md at
# the top level directory of deal.II.
#
# ---------------------------------------------------------------------

import math
import unittest
from PyDealII.Debug import *

class TestManifoldWrapperShell(unittest.TestCase):

    def setUp(self):
        self.triangulation = Triangulation('2D')
        p_center = Point([0, 0])
        self.triangulation.generate_hyper_shell(center = p_center, inner_radius = 0.5, outer_radius = 1., n_cells = 0, colorize = True)

        self.manifold = Manifold(dim = 2, spacedim = 2)
        self.manifold.create_polar(p_center)

    def test_manifold(self):
        self.triangulation.set_manifold(0, self.manifold)
        for cell in self.triangulation.active_cells():
            cell.manifold_id = 0

        self.triangulation.refine_global(3)

        circumference = 0
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary() and face.boundary_id == 1:
                        circumference += face.measure()

        self.assertTrue(abs(circumference - 2*math.pi)/(2*math.pi) < 1e-2)

class TestManifoldWrapperBall(unittest.TestCase):

    def setUp(self):
        self.triangulation = Triangulation('3D')
        p_center = Point([0., 0., 0.])
        self.triangulation.generate_hyper_ball(center = p_center, radius = 1.)

        self.manifold = Manifold(dim = 3, spacedim = 3)
        self.manifold.create_spherical(p_center)

    def test_manifold(self):
        self.triangulation.reset_manifold(number = 0)
        self.triangulation.set_manifold(number = 0, manifold = self.manifold)
        for cell in self.triangulation.active_cells():
            if cell.at_boundary():
                cell.manifold_id = 0

        self.triangulation.refine_global(3)

        volume = 0
        for cell in self.triangulation.active_cells():
            volume += cell.measure()

        self.assertTrue(abs(volume - 4./3. * math.pi) / (4./3.*math.pi) < 2e-2)

if __name__ == '__main__':
    unittest.main()
