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


class TestTriangulationWrapper(unittest.TestCase):

    def setUp(self):
        self.dim = [['2D', '2D'], ['2D', '3D'], ['3D', '3D']]

    def build_hyper_cube_triangulation(self, dim):
        triangulation_1 = Triangulation(dim[0], dim[1])
        triangulation_1.generate_hyper_cube()
        return triangulation_1

    def build_hyper_rectangle_triangulation(self, dim):
        triangulation_2 = Triangulation(dim[0], dim[1])
        if (dim[0] == '2D'):
            point_1 = Point([0., 0.])
            point_2 = Point([1., 1.])
        else:
            point_1 = Point([0., 0., 0.])
            point_2 = Point([1., 1., 1.])
        triangulation_2.generate_hyper_rectangle(point_1, point_2)
        return triangulation_2

    def test_hyper_cube(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            n_cells = triangulation_1.n_active_cells()
            self.assertEqual(n_cells, 1)

    def test_simplex(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if (dim[0] == '2D'):
                point_1 = Point([0., 0.])
                point_2 = Point([1., 0.])
                point_3 = Point([1., 1.])
                vertices = [point_1, point_2, point_3]
            else:
                point_1 = Point([0., 0., 0.])
                point_2 = Point([1., 0., 0.])
                point_3 = Point([1., 1., 0.])
                point_4 = Point([1., 1., 1.])
                vertices = [point_1, point_2, point_3, point_4]
            triangulation.generate_simplex(vertices)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, len(vertices))

    def test_subdivided_hyper_cube(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0], dim[1])
            repetitions = 2
            triangulation.generate_subdivided_hyper_cube(repetitions)
            n_cells = triangulation.n_active_cells()
            if (dim[0] == '2D'):
                self.assertEqual(n_cells, 4)
            else:
                self.assertEqual(n_cells, 8)

    def test_hyper_rectangle(self):
        for dim in self.dim:
            triangulation_2 = self.build_hyper_rectangle_triangulation(dim)
            n_cells = triangulation_2.n_active_cells()
            self.assertEqual(n_cells, 1)

    def test_subdivided_hyper_rectangle(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0], dim[1])
            if (dim[0] == '2D'):
                repetitions = [6, 4]
                point_1 = Point([0., 0.])
                point_2 = Point([1., 1.])
            else:
                repetitions = [2, 3, 4]
                point_1 = Point([0., 0., 0.])
                point_2 = Point([1., 1., 1.])
            triangulation.generate_subdivided_hyper_rectangle(repetitions,
                                                              point_1,
                                                              point_2)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 24)

    def test_hyper_ball(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if (dim[0] == '2D'):
                center = Point([0., 0.])
                n_cells_ref = 5
            else:
                center = Point([0., 0., 0.])
                n_cells_ref = 7
            triangulation.generate_hyper_ball(center)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_shift_and_merge(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_2 = self.build_hyper_rectangle_triangulation(dim)
            if (dim[1] == '2D'):
                triangulation_2.shift([1., 0.])
            else:
                triangulation_2.shift([1., 0., 0.])
            triangulation = Triangulation(dim[0], dim[1])
            triangulation.merge_triangulations(triangulation_1,
                                               triangulation_2)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 2)

    def test_adaptive_refinement(self):
        for dim in self.dim:
            triangulation = self.build_hyper_cube_triangulation(dim)
            triangulation.refine_global(1)
            for cell in triangulation.active_cells():
                cell.refine_flag = 'isotropic'
                break
            triangulation.execute_coarsening_and_refinement()
            n_cells = triangulation.n_active_cells()
            if (dim[0] == '2D'):
                self.assertEqual(n_cells, 7)
            else:
                self.assertEqual(n_cells, 15)

    def test_refine_global(self):
        for dim in self.dim:
            triangulation = self.build_hyper_cube_triangulation(dim)
            triangulation.refine_global(2)
            n_cells = triangulation.n_active_cells()
            if (dim[0] == '2D'):
                self.assertEqual(n_cells, 16)
            else:
                self.assertEqual(n_cells, 64)

    def test_save_load(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_1.refine_global(1)
            triangulation_1.save('mesh.output')
            triangulation_2 = Triangulation(dim[0], dim[1])
            triangulation_2.load('mesh.output')
            n_cells = triangulation_2.n_active_cells()
            if (dim[0] == '2D'):
                self.assertEqual(n_cells, 4)
            else:
                self.assertEqual(n_cells, 8)


if __name__ == '__main__':
    unittest.main()
