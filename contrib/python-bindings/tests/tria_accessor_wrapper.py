## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2024 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

import math
import unittest

try:
    from PyDealII.Debug import *
except ImportError:
    from PyDealII.Release import *


class TestTriaAccessorWrapper2D(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("2D")
        self.triangulation.generate_hyper_cube()
        self.triangulation.refine_global(1)

    def test_at_boundary(self):
        # On a once-refined unit square, boundary faces are on the outer edges.
        # Count boundary faces and interior faces separately.
        n_boundary = 0
        n_interior = 0
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    n_boundary += 1
                else:
                    n_interior += 1
        # 4 cells, each with 4 faces = 16 face-slots.
        # 8 boundary half-edges + 4 interior edges counted twice (once per cell) = 8+8=16.
        self.assertEqual(n_boundary, 8)
        self.assertEqual(n_interior, 8)

    def test_boundary_id(self):
        # By default all boundary faces have boundary_id 0. Verify, then change.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    self.assertEqual(face.boundary_id, 0)

        # Set each boundary face to a distinct id.
        bid = 1
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    face.boundary_id = bid
                    bid += 1

        # Verify that no boundary face still has id 0.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    self.assertNotEqual(face.boundary_id, 0)

    def test_manifold_id(self):
        mid = 5
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                face.manifold_id = mid
                self.assertEqual(face.manifold_id, mid)
                mid += 1

    def test_set_all_boundary_ids(self):
        # set_all_boundary_ids propagates the id to all sub-objects of the face.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    face.set_all_boundary_ids(42)
                    self.assertEqual(face.boundary_id, 42)

    def test_measure(self):
        # After one global refinement the unit square contains 4 cells of side
        # length 0.5, so each face (edge) has measure 0.5.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                self.assertAlmostEqual(face.measure(), 0.5)

        # Total boundary length must equal the perimeter of the unit square.
        total_boundary_length = 0.0
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    total_boundary_length += face.measure()
        self.assertAlmostEqual(total_boundary_length, 4.0)

    def test_barycenter(self):
        # Verify that face barycenters lie on the unit-square boundary or interior.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                bc = face.barycenter()
                self.assertTrue(0.0 <= bc.x <= 1.0)
                self.assertTrue(0.0 <= bc.y <= 1.0)

    def test_center_matches_barycenter_for_flat_faces(self):
        # For flat (straight) faces center() and barycenter() must agree.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                bc = face.barycenter()
                c = face.center()
                self.assertAlmostEqual(bc.x, c.x)
                self.assertAlmostEqual(bc.y, c.y)

    def test_get_vertex(self):
        # Each 2-D face (edge) has exactly 2 vertices.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                v0 = face.get_vertex(0)
                v1 = face.get_vertex(1)
                # Both vertices must be inside (or on the boundary of) [0, 1]^2.
                for v in (v0, v1):
                    self.assertTrue(0.0 <= v.x <= 1.0)
                    self.assertTrue(0.0 <= v.y <= 1.0)
                # The two endpoints must be distinct.
                self.assertFalse(v0.x == v1.x and v0.y == v1.y)

    def test_set_vertex(self):
        # Move one vertex of the very first face and verify it changed.
        for cell in self.triangulation.active_cells():
            face = cell.faces()[0]
            original = face.get_vertex(0)
            new_point = Point([original.x + 0.05, original.y + 0.05])
            face.set_vertex(0, new_point)
            moved = face.get_vertex(0)
            self.assertAlmostEqual(moved.x, new_point.x)
            self.assertAlmostEqual(moved.y, new_point.y)
            break


class TestTriaAccessorWrapper3D(unittest.TestCase):
    def setUp(self):
        self.triangulation = Triangulation("3D")
        self.triangulation.generate_hyper_cube()

    def test_at_boundary(self):
        # An unrefined unit cube has 1 cell with 6 boundary faces.
        n_boundary = 0
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    n_boundary += 1
        self.assertEqual(n_boundary, 6)

    def test_boundary_id(self):
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                if face.at_boundary():
                    self.assertEqual(face.boundary_id, 0)
                    face.boundary_id = 7
                    self.assertEqual(face.boundary_id, 7)

    def test_manifold_id(self):
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                face.manifold_id = 3
                self.assertEqual(face.manifold_id, 3)

    def test_measure(self):
        # Each face of the unrefined unit cube is a unit square with area 1.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                self.assertAlmostEqual(face.measure(), 1.0)

    def test_measure_after_refinement(self):
        # After one global refinement each face has area 0.25.
        self.triangulation.refine_global(1)
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                self.assertAlmostEqual(face.measure(), 0.25)

    def test_barycenter(self):
        for cell in self.triangulation.active_cells():
            bc = cell.barycenter()
            self.assertTrue(0.0 <= bc.x <= 1.0)
            self.assertTrue(0.0 <= bc.y <= 1.0)
            self.assertTrue(0.0 <= bc.z <= 1.0)

    def test_get_vertex(self):
        # Each 3-D face (quad) has exactly 4 vertices.
        for cell in self.triangulation.active_cells():
            for face in cell.faces():
                vertices = [face.get_vertex(i) for i in range(4)]
                for v in vertices:
                    self.assertTrue(0.0 <= v.x <= 1.0)
                    self.assertTrue(0.0 <= v.y <= 1.0)
                    self.assertTrue(0.0 <= v.z <= 1.0)


if __name__ == "__main__":
    unittest.main()
