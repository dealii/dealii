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

class TestCellAccessorWrapper(unittest.TestCase):

    def setUp(self):
        self.triangulation = Triangulation('2D')
        self.triangulation.generate_hyper_cube()
        self.triangulation.refine_global(1)

    def test_faces(self):
        for cell in self.triangulation.active_cells():
            faces = cell.faces()
            self.assertEqual(len(faces), 4)

    def test_at_boundary(self):
        for cell in self.triangulation.active_cells():
            self.assertEqual(cell.at_boundary(), cell.has_boundary_lines()) 

    def test_faces(self):
        n_neighbors = 0
        for cell in self.triangulation.active_cells():
            faces = cell.faces()
            for i in range(len(faces)):
                if not faces[i].at_boundary():
                    neighbor = cell.neighbor(i)
                    n_neighbors += 1
        self.assertEqual(n_neighbors, 8)

    def test_material_id(self):
        material_id = 0
        for cell in self.triangulation.active_cells():
            cell.material_id = material_id
            material_id += 1
        material_id = 0
        for cell in self.triangulation.active_cells():
            self.assertEqual(cell.material_id, material_id)
            material_id += 1

    def test_manifold_id(self):
        manifold_id = 0
        for cell in self.triangulation.active_cells():
            cell.manifold_id = manifold_id
            manifold_id += 1
        manifold_id = 0
        for cell in self.triangulation.active_cells():
            self.assertEqual(cell.manifold_id, manifold_id)
            manifold_id += 1

    def test_refine_flag(self):
        index = 0
        refine_flags = ['no_refinement', 'cut_x', 'cut_y', 'cut_xy']
        for cell in self.triangulation.active_cells():
            cell.refine_flag = refine_flags[index]
            index += 1
        index = 0
        for cell in self.triangulation.active_cells():
            self.assertEqual(cell.refine_flag, refine_flags[index])
            index += 1

    def test_coarsen_flag(self):
        coarsen_flag = True
        for cell in self.triangulation.active_cells():
            cell.coarsen_flag = coarsen_flag
            coarsen_flag = not coarsen_flag
        coarsen_flag = True
        for cell in self.triangulation.active_cells():
            self.assertEqual(cell.coarsen_flag, coarsen_flag)
            coarsen_flag = not coarsen_flag

    def test_barycenter(self):
        centers = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
        index = 0
        for cell in self.triangulation.active_cells():
            barycenter = cell.barycenter()
            self.assertEqual(barycenter.x, centers[index][0])
            self.assertEqual(barycenter.y, centers[index][1])
            index += 1

    def test_move_vertex(self):
        point = Point([0.6, 0.6])
        for cell in self.triangulation.active_cells():
            cell.set_vertex(3, point)
            vertex = cell.get_vertex(3)
            break
        vertices = [3, 2, 1, 0]
        index = 0
        for cell in self.triangulation.active_cells():
            vertex = cell.get_vertex(vertices[index])
            self.assertEqual(vertex.x, point.x)
            self.assertEqual(vertex.y, point.y)
            index += 1


if __name__ == '__main__':
    unittest.main()
