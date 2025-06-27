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

import math
import os
import unittest

try:
    from PyDealII.Debug import *
except ImportError:
    from PyDealII.Release import *


class TestManifoldWrapperShell(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("2D")
        p_center = Point([0, 0])
        self.triangulation.generate_hyper_shell(
            center=p_center,
            inner_radius=0.5,
            outer_radius=1.0,
            n_cells=0,
            colorize=True,
        )

        self.manifold = Manifold(dim=2, spacedim=2)
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

        self.assertTrue(abs(circumference - 2 * math.pi) / (2 * math.pi) < 1e-2)


class TestManifoldWrapperBall(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("3D")
        p_center = Point([0.0, 0.0, 0.0])
        self.triangulation.generate_hyper_ball(center=p_center, radius=1.0)

        self.manifold = Manifold(dim=3, spacedim=3)
        self.manifold.create_spherical(p_center)

    def test_manifold(self):
        self.triangulation.reset_manifold(number=0)
        self.triangulation.set_manifold(number=0, manifold=self.manifold)
        for cell in self.triangulation.active_cells():
            if cell.at_boundary():
                cell.manifold_id = 0

        self.triangulation.refine_global(3)

        volume = 0
        for cell in self.triangulation.active_cells():
            volume += cell.measure()

        self.assertTrue(
            abs(volume - 4.0 / 3.0 * math.pi) / (4.0 / 3.0 * math.pi) < 2e-2
        )


class TestManifoldWrapperFunction(unittest.TestCase):
    def setUp(self):
        self.manifold_1 = Manifold(dim=2, spacedim=2)
        self.manifold_1.create_function_string("x^2;y^2", "sqrt(x);sqrt(y)")

        self.manifold_2 = Manifold(dim=2, spacedim=2)
        self.manifold_2.create_function(
            lambda p: [p[0] ** 2.0, p[1] ** 2.0],
            lambda p: [math.sqrt(p[0]), math.sqrt(p[1])],
        )

        self.tria_reference = Triangulation("2D")
        test_directory = os.environ.get("DEAL_II_PYTHON_TESTPATH")
        self.tria_reference.read(test_directory + "/manifold_wrapper.vtk", "vtk")

    def test_manifold_str(self):
        self.triangulation = Triangulation("2D")
        p_center = Point([0.0, 0.0, 0.0])
        self.triangulation.generate_hyper_cube()
        self.triangulation.reset_manifold(number=0)
        self.triangulation.set_manifold(number=0, manifold=self.manifold_1)
        for cell in self.triangulation.active_cells():
            cell.set_all_manifold_ids(0)

        self.triangulation.refine_global(2)

        for cell_ref, cell in zip(
            self.tria_reference.active_cells(), self.triangulation.active_cells()
        ):
            self.assertTrue(abs(cell_ref.measure() - cell.measure()) < 1e-8)

    def test_manifold_lambda(self):
        self.triangulation = Triangulation("2D")
        p_center = Point([0.0, 0.0, 0.0])
        self.triangulation.generate_hyper_cube()
        self.triangulation.reset_manifold(number=0)
        self.triangulation.set_manifold(number=0, manifold=self.manifold_2)
        for cell in self.triangulation.active_cells():
            cell.set_all_manifold_ids(0)

        self.triangulation.refine_global(2)

        for cell_ref, cell in zip(
            self.tria_reference.active_cells(), self.triangulation.active_cells()
        ):
            self.assertTrue(abs(cell_ref.measure() - cell.measure()) < 1e-8)


if __name__ == "__main__":
    unittest.main()
