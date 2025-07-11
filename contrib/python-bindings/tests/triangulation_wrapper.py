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

try:
    from PyDealII.Debug import *
except ImportError:
    from PyDealII.Release import *


class TestTriangulationWrapper(unittest.TestCase):
    def setUp(self):
        self.dim = [["2D", "2D"], ["2D", "3D"], ["3D", "3D"]]
        self.restricted_dim = [["2D", "2D"], ["3D", "3D"]]

    def build_hyper_cube_triangulation(self, dim):
        triangulation_1 = Triangulation(dim[0], dim[1])
        triangulation_1.generate_hyper_cube()
        return triangulation_1

    def build_hyper_rectangle_triangulation(self, dim):
        triangulation_2 = Triangulation(dim[0], dim[1])
        if dim[0] == "2D":
            point_1 = Point([0.0, 0.0])
            point_2 = Point([1.0, 1.0])
        else:
            point_1 = Point([0.0, 0.0, 0.0])
            point_2 = Point([1.0, 1.0, 1.0])
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
            if dim[0] == "2D":
                point_1 = Point([0.0, 0.0])
                point_2 = Point([1.0, 0.0])
                point_3 = Point([1.0, 1.0])
                vertices = [point_1, point_2, point_3]
            else:
                point_1 = Point([0.0, 0.0, 0.0])
                point_2 = Point([1.0, 0.0, 0.0])
                point_3 = Point([1.0, 1.0, 0.0])
                point_4 = Point([1.0, 1.0, 1.0])
                vertices = [point_1, point_2, point_3, point_4]
            triangulation.generate_simplex(vertices)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, len(vertices))

    def test_create_triangulation(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
                cell_vertices = [[0, 1, 2, 3]]
            else:
                vertices = [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [0.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0],
                ]
                cell_vertices = [[0, 1, 2, 3, 4, 5, 6, 7]]
            triangulation.create_triangulation(vertices, cell_vertices)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, len(cell_vertices))

    def test_subdivided_hyper_cube(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0], dim[1])
            repetitions = 2
            triangulation.generate_subdivided_hyper_cube(repetitions)
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
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
            if dim[0] == "2D":
                repetitions = [6, 4]
                point_1 = Point([0.0, 0.0])
                point_2 = Point([1.0, 1.0])
            else:
                repetitions = [2, 3, 4]
                point_1 = Point([0.0, 0.0, 0.0])
                point_2 = Point([1.0, 1.0, 1.0])
            triangulation.generate_subdivided_hyper_rectangle(
                repetitions, point_1, point_2
            )
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 24)

    def test_subdivided_steps_hyper_rectangle(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0], dim[1])
            if dim[0] == "2D":
                step_sizes = [[0.6, 0.4], [0.4, 0.6]]
                point_1 = Point([0.0, 0.0])
                point_2 = Point([1.0, 1.0])
                triangulation.generate_subdivided_steps_hyper_rectangle(
                    step_sizes, point_1, point_2
                )
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 4)
            else:
                step_sizes = [[0.6, 0.4], [0.4, 0.6], [0.5, 0.5]]
                point_1 = Point([0.0, 0.0, 0.0])
                point_2 = Point([1.0, 1.0, 1.0])
                triangulation.generate_subdivided_steps_hyper_rectangle(
                    step_sizes, point_1, point_2
                )
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 8)

    def test_subdivided_material_hyper_rectangle(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0], dim[1])
            if dim[0] == "2D":
                spacing = [[0.3, 0.3, 0.4], [0.4, 0.3, 0.3]]
                point = Point([1.0, 1.0])
                material_ids = [[1, 2, 3], [0, -1, 1], [2, -1, 3]]
                triangulation.generate_subdivided_material_hyper_rectangle(
                    spacing, point, material_ids
                )
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 7)
            else:
                spacing = [[0.3, 0.3, 0.4], [0.4, 0.3, 0.3], [0.2, 0.3, 0.5]]
                point = Point([1.0, 1.0, 1.0])
                material_ids = [
                    [[1, 2, 3], [0, -1, 1], [2, -1, 3]],
                    [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                    [[2, -1, -1], [-1, -1, -1], [-1, -1, -1]],
                ]
                triangulation.generate_subdivided_material_hyper_rectangle(
                    spacing, point, material_ids
                )
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 17)

    def test_hyper_cube_with_cylindrical_hole(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0])
            triangulation.generate_hyper_cube_with_cylindrical_hole(
                inner_radius=0.25,
                outer_radius=0.5,
                L=0.5,
                repetitions=1,
                colorize=False,
            )
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 8)

    def test_generate_cheese(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0], dim[1])
            if dim[0] == "2D":
                holes = [2, 3]
                triangulation.generate_cheese(holes)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 29)
            else:
                holes = [2, 3, 2]
                triangulation.generate_cheese(holes)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 175)

    def test_generate_plate_with_hole(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            triangulation.generate_plate_with_a_hole()
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 28)

    def test_generate_channel_with_cylinder(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            triangulation.generate_channel_with_cylinder()
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 108)
            else:
                self.assertEqual(n_cells, 432)

    def test_generate_general_cell(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0], dim[1])
            if dim[0] == "2D":
                point_0 = Point([0.0, 0.0])
                point_1 = Point([1.0, 0.0])
                point_2 = Point([0.0, 1.0])
                point_3 = Point([1.0, 1.0])
                points = [point_0, point_1, point_2, point_3]
                triangulation.generate_general_cell(points)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 4)
            else:
                point_0 = Point([0.0, 0.0, 0.0])
                point_1 = Point([1.0, 0.0, 0.0])
                point_2 = Point([0.0, 1.0, 0.0])
                point_3 = Point([1.0, 1.0, 0.0])
                point_4 = Point([0.0, 0.0, 1.0])
                point_5 = Point([1.0, 0.0, 1.0])
                point_6 = Point([0.0, 1.0, 1.0])
                point_7 = Point([1.0, 1.0, 1.0])
                points = [
                    point_0,
                    point_1,
                    point_2,
                    point_3,
                    point_4,
                    point_5,
                    point_6,
                    point_7,
                ]
                triangulation.generate_general_cell(points)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 8)

    def test_generate_parallelogram(self):
        triangulation = Triangulation("2D")
        corner_0 = Point([1.0, 0.0])
        corner_1 = Point([1.0, 1.0])
        corners = [corner_0, corner_1]
        triangulation.generate_parallelogram(corners)
        triangulation.refine_global(1)
        n_cells = triangulation.n_active_cells()
        self.assertEqual(n_cells, 4)

    def test_generate_parallelepid(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                corner_0 = Point([1.0, 0.0])
                corner_1 = Point([1.0, 1.0])
                corners = [corner_0, corner_1]
                triangulation.generate_parallelepiped(corners)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 4)
            else:
                corner_0 = Point([1.0, 0.0, 0.0])
                corner_1 = Point([0.0, 1.0, 0.0])
                corner_2 = Point([0.0, 0.0, 1.0])
                corners = [corner_0, corner_1, corner_2]
                triangulation.generate_parallelepiped(corners)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 8)

    def test_generate_fixed_subdivided_parallelepiped(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                corner_0 = Point([1.0, 0.0])
                corner_1 = Point([1.0, 1.0])
                corners = [corner_0, corner_1]
                triangulation.generate_fixed_subdivided_parallelepiped(2, corners)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 16)
            else:
                corner_0 = Point([1.0, 0.0, 0.0])
                corner_1 = Point([0.0, 1.0, 0.0])
                corner_2 = Point([0.0, 0.0, 1.0])
                corners = [corner_0, corner_1, corner_2]
                triangulation.generate_fixed_subdivided_parallelepiped(2, corners)
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 64)

    def test_generate_varying_subdivided_parallelepiped(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                corner_0 = Point([1.0, 0.0])
                corner_1 = Point([1.0, 1.0])
                corners = [corner_0, corner_1]
                subdivisions = [2, 3]
                triangulation.generate_varying_subdivided_parallelepiped(
                    subdivisions, corners
                )
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 24)
            else:
                corner_0 = Point([1.0, 0.0, 0.0])
                corner_1 = Point([0.0, 1.0, 0.0])
                corner_2 = Point([0.0, 0.0, 1.0])
                corners = [corner_0, corner_1, corner_2]
                subdivisions = [2, 3, 1]
                triangulation.generate_varying_subdivided_parallelepiped(
                    subdivisions, corners
                )
                triangulation.refine_global(1)
                n_cells = triangulation.n_active_cells()
                self.assertEqual(n_cells, 48)

    def test_generate_enclosed_hyper_cube(self):
        for dim in self.restricted_dim:
            triangulation = Triangulation(dim[0])
            left = 1.0
            right = 3.0
            thickness = 4.0
            triangulation.generate_enclosed_hyper_cube(left, right, thickness)
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 9)
            else:
                self.assertEqual(n_cells, 27)

    def test_hyper_ball(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                center = Point([0.0, 0.0])
                n_cells_ref = 5
            else:
                center = Point([0.0, 0.0, 0.0])
                n_cells_ref = 7
            triangulation.generate_hyper_ball(center)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_hyper_ball_balanced(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            triangulation.generate_hyper_ball_balanced()
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 12)
            else:
                self.assertEqual(n_cells, 32)

    def test_hyper_sphere(self):
        triangulation = Triangulation("2D", "3D")
        center = Point([0, 0])
        radius = 1.0
        triangulation.generate_hyper_sphere(center, radius)
        n_cells = triangulation.n_active_cells()
        self.assertEqual(n_cells, 6)

    def test_quarter_hyper_ball(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                center = Point([0.0, 0.0])
                n_cells_ref = 3
            else:
                center = Point([0.0, 0.0, 0.0])
                n_cells_ref = 4
            triangulation.generate_quarter_hyper_ball(center)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_half_hyper_ball(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                center = Point([0.0, 0.0])
                n_cells_ref = 4
            else:
                center = Point([0.0, 0.0, 0.0])
                n_cells_ref = 6
            triangulation.generate_half_hyper_ball(center)
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_cylinder(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                n_cells_ref = 1
            else:
                n_cells_ref = 10
            triangulation.generate_cylinder()
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_subdivided_cylinder(self):
        triangulation = Triangulation("3D")
        n_cells_ref = 5
        triangulation.generate_subdivided_cylinder(2)
        n_cells = triangulation.n_active_cells()
        self.assertEqual(n_cells, n_cells_ref)

    def test_truncated_cone(self):
        for dim in self.dim:
            triangulation = Triangulation(dim[0])
            if dim[0] == "2D":
                n_cells_ref = 1
            else:
                n_cells_ref = 5
            triangulation.generate_truncated_cone()
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, n_cells_ref)

    def test_shift_and_merge(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_2 = self.build_hyper_rectangle_triangulation(dim)
            if dim[1] == "2D":
                triangulation_2.shift([1.0, 0.0])
            else:
                triangulation_2.shift([1.0, 0.0, 0.0])
            triangulation = Triangulation(dim[0], dim[1])
            triangulation.merge_triangulations([triangulation_1, triangulation_2])
            n_cells = triangulation.n_active_cells()
            self.assertEqual(n_cells, 2)

    def test_replicate(self):
        for dim in self.restricted_dim:
            triangulation_in = self.build_hyper_cube_triangulation(dim)
            triangulation_out = Triangulation(dim[0])
            if dim[0] == "2D":
                triangulation_out.replicate_triangulation(triangulation_in, [3, 2])
            else:
                triangulation_out.replicate_triangulation(triangulation_in, [3, 2, 1])
            n_cells = triangulation_out.n_active_cells()
            self.assertEqual(n_cells, 6)

    def test_flatten(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_1.refine_global(2)
            triangulation_2 = Triangulation(dim[0])
            triangulation_1.flatten_triangulation(triangulation_2)
            n_cells = triangulation_2.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 16)
            else:
                self.assertEqual(n_cells, 64)

    def test_adaptive_refinement(self):
        for dim in self.dim:
            triangulation = self.build_hyper_cube_triangulation(dim)
            triangulation.refine_global(1)
            for cell in triangulation.active_cells():
                cell.refine_flag = "isotropic"
                break
            triangulation.execute_coarsening_and_refinement()
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 7)
            else:
                self.assertEqual(n_cells, 15)

    def test_refine_global(self):
        for dim in self.dim:
            triangulation = self.build_hyper_cube_triangulation(dim)
            triangulation.refine_global(2)
            n_cells = triangulation.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 16)
            else:
                self.assertEqual(n_cells, 64)

    def test_transform(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_1.refine_global(1)
            triangulation_2 = self.build_hyper_cube_triangulation(dim)
            triangulation_2.refine_global(1)

            triangulation_1.transform(lambda p: [v + 1.0 for v in p])

            if dim[1] == "3D":
                offset = Point([1.0, 1.0, 1.0])
            else:
                offset = Point([1.0, 1.0])

            for (cell_1, cell_2) in zip(
                triangulation_1.active_cells(), triangulation_2.active_cells()
            ):
                self.assertTrue(
                    cell_1.center().distance(cell_2.center() + offset) < 1e-8
                )

    def test_find_active_cell_around_point(self):
        for dim in self.dim:
            triangulation = self.build_hyper_cube_triangulation(dim)
            triangulation.refine_global(2)

            for cell in triangulation.active_cells():
                cell_ret = triangulation.find_active_cell_around_point(cell.center())
                self.assertTrue(cell.center().distance(cell_ret.center()) < 1e-8)

    def test_simplex(self):
        for dim in self.dim:
            triangulation_hex = self.build_hyper_cube_triangulation(dim)
            triangulation_simplex = Triangulation(dim[0], dim[1])
            triangulation_hex.convert_hypercube_to_simplex_mesh(triangulation_simplex)

            if dim[0] == "3D":
                self.assertTrue(triangulation_simplex.n_active_cells() == 24)
            else:
                self.assertTrue(triangulation_simplex.n_active_cells() == 8)

    def test_save_load(self):
        for dim in self.dim:
            triangulation_1 = self.build_hyper_cube_triangulation(dim)
            triangulation_1.refine_global(1)
            triangulation_1.save("mesh.output")
            triangulation_2 = Triangulation(dim[0], dim[1])
            triangulation_2.load("mesh.output")
            n_cells = triangulation_2.n_active_cells()
            if dim[0] == "2D":
                self.assertEqual(n_cells, 4)
            else:
                self.assertEqual(n_cells, 8)

    def test_mesh_smoothing(self):
        tria = Triangulation("2D")
        self.assertEqual(tria.get_mesh_smoothing(), MeshSmoothing.none)
        tria.set_mesh_smoothing(MeshSmoothing.maximum_smoothing)
        self.assertEqual(tria.get_mesh_smoothing(), MeshSmoothing.maximum_smoothing)

        tria = Triangulation("2D", MeshSmoothing.limit_level_difference_at_vertices)
        self.assertEqual(
            tria.get_mesh_smoothing(), MeshSmoothing.limit_level_difference_at_vertices
        )

        tria = Triangulation("2D", "3D", MeshSmoothing.none, False)
        tria.set_mesh_smoothing(MeshSmoothing.limit_level_difference_at_vertices)
        self.assertEqual(
            tria.get_mesh_smoothing(), MeshSmoothing.limit_level_difference_at_vertices
        )

        tria = Triangulation(
            "3D",
            MeshSmoothing.limit_level_difference_at_vertices
            | MeshSmoothing.do_not_produce_unrefined_islands,
        )
        self.assertEqual(
            tria.get_mesh_smoothing(),
            MeshSmoothing.limit_level_difference_at_vertices
            | MeshSmoothing.do_not_produce_unrefined_islands,
        )


if __name__ == "__main__":
    unittest.main()
