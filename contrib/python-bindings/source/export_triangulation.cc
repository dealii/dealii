// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <boost/python.hpp>

#include <cell_accessor_wrapper.h>
#include <triangulation_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  // Macro to enable default arguments
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_cube_overloads,
                                         generate_hyper_cube,
                                         0,
                                         3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_subdivided_hyper_cube_overloads,
    generate_subdivided_hyper_cube,
    1,
    3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_rectangle_overloads,
                                         generate_hyper_rectangle,
                                         2,
                                         3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_subdivided_hyper_rectangle_overloads,
    generate_subdivided_hyper_rectangle,
    3,
    4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_subdivided_steps_hyper_rectangle_overloads,
    generate_subdivided_steps_hyper_rectangle,
    3,
    4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_subdivided_material_hyper_rectangle_overloads,
    generate_subdivided_material_hyper_rectangle,
    3,
    4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_general_cell_overloads,
                                         generate_general_cell,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_parallelogram_overloads,
                                         generate_parallelogram,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_parallelepiped_overloads,
                                         generate_parallelepiped,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_fixed_subdivided_parallelepiped_overloads,
    generate_fixed_subdivided_parallelepiped,
    2,
    3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_varying_subdivided_parallelepiped_overloads,
    generate_varying_subdivided_parallelepiped,
    2,
    3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_enclosed_hyper_cube_overloads,
                                         generate_enclosed_hyper_cube,
                                         0,
                                         4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_hyper_cube_with_cylindrical_hole_overloads,
    generate_hyper_cube_with_cylindrical_hole,
    0,
    5)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_plate_with_a_hole_overloads,
                                         generate_plate_with_a_hole,
                                         0,
                                         12)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    generate_channel_with_cylinder_overloads,
    generate_channel_with_cylinder,
    0,
    4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_ball_overloads,
                                         generate_hyper_ball,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_ball_balanced_overloads,
                                         generate_hyper_ball_balanced,
                                         0,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_sphere_overloads,
                                         generate_hyper_sphere,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_quarter_hyper_ball_overloads,
                                         generate_quarter_hyper_ball,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_half_hyper_ball_overloads,
                                         generate_half_hyper_ball,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_cylinder_overloads,
                                         generate_cylinder,
                                         0,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_subdivided_cylinder_overloads,
                                         generate_subdivided_cylinder,
                                         1,
                                         3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_truncated_cone_overloads,
                                         generate_truncated_cone,
                                         0,
                                         3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_shell_overloads,
                                         generate_hyper_shell,
                                         3,
                                         5)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(distort_random_overloads,
                                         distort_random,
                                         1,
                                         2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    find_active_cell_around_point_overloads,
    find_active_cell_around_point,
    1,
    2)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(merge_triangulations_overloads,
                                         merge_triangulations,
                                         1,
                                         3)

  const char n_active_cells_docstring[] =
    "Return the number of active cells.                                     \n";



  const char n_cells_docstring[] =
    "Return the number of cells.                                            \n";



  const char dim_docstring[] =
    "Return the dimension of the Triangulation                              \n";



  const char spacedim_docstring[] =
    "Return the space dimension of the Triangulation.                       \n";



  const char create_triangulation_docstring[] =
    "Given a list of points and how vertices connect to cells, create a     \n"
    "mesh.                                                                  \n";



  const char generate_hyper_cube_docstring[] =
    "Generate a hyper_cube (square in 2D and cube in 3D)                    \n"
    "with exactly one cell.                                                 \n";



  const char generate_simplex_docstring[] =
    "Generate a simplex with (dim+1) vertices and mesh cells.               \n";



  const char generate_subdivided_hyper_cube_docstring[] =
    "Same as hyper_cube but not only one cell is created but each           \n"
    "coordinate direction is subdivided in repetitions cells                \n";



  const char generate_hyper_rectangle_docstring[] =
    "Generate a coordinate-parallel brick from the two diagonally opposite  \n"
    "corners points p1 and p2.                                              \n";



  const char generate_subdivided_hyper_rectangle_docstring[] =
    "Generate a coordinate-parallel brick from the two diagonally opposite  \n"
    "corners point p1 and p2. In direction i, repetitions[i] cells are      \n"
    "created.                                                               \n";



  const char generate_subdivided_steps_hyper_rectangle_docstring[] =
    "Like the previous function. However, here the first argument does not  \n"
    "denote the number of subdivisions in each coordinate direction, but a  \n"
    "sequence of step sizes for each coordinate direction. This function is \n"
    "therefore the right one to generate graded meshes where cells are      \n"
    "concentrated in certains areas, rather than a uniformly subdidived mesh\n"
    "as the previous function generates.                                    \n";



  const char generate_subdivided_material_hyper_rectangle_docstring[] =
    "Like the previous function, but with the following twist: the          \n"
    "material_id argument is a dim-dimensional array that, for each cell,   \n"
    "indicates which material_id should be set. In addition, and this is the\n"
    "major new functionality, if the material_id of a cell is (-1), then    \n"
    "that cell is deleted from the triangulation, i.e. the domain will have \n"
    "a void there.                                                          \n";



  const char generate_cheese_docstring[] =
    "Rectangular domain with rectangular pattern of holes. The domain itself\n"
    "is rectangular, very much as if it had been generated by               \n"
    "subdivided_hyper_rectangle(). The argument holes specifies how many    \n"
    "square holes the domain should have in each coordinate direction. The  \n"
    "total number of mesh cells in that direction is then this number plus  \n"
    "one. The number of holes in one direction must be at least one.        \n";



  const char generate_general_cell_docstring[] =
    "A general quadrilateral in 2d or a general hexahedron in 3d. It is the \n"
    "responsibility of the user to provide the vertices in the right order  \n"
    "(see the documentation of the GeometryInfo class) because the vertices \n"
    "are stored in the same order as they are given. It is also important to\n"
    "make that the volume of the cell is positive. If the argument          \n"
    "colorize is false, all boundary indicators are set to zero (not        \n"
    "colorized) for 2d and 3d. If it is true, the boundary is colorized as  \n"
    "in hyper_rectangle(). In 1d, the indicators are always colorized.      \n";



  const char generate_parallelogram_docstring[] =
    "A parallelogram. The first corner point is the origin. The dim         \n"
    "adjacent points are the ones given in the second argument and the      \n"
    "fourth point will be the sum of these two vectors. Colorizing is done  \n"
    "in the same way as in hyper_rectangle().                               \n"
    "Note: This function is implemented in 2d only.                         \n";



  const char generate_parallelepiped_docstring[] =
    "A parallelepiped. The first corner point is the origin. The dim        \n"
    "adjacent points are vectors describing the edges of the parallelepiped \n"
    "with respect to the origin. Additional points are sums of these dim    \n"
    "vectors. Colorizing is done according to hyper_rectangle().            \n"
    "Note: This function silently reorders the vertices on the cells to     \n"
    "lexicographic ordering (see GridTools::consistently_order_cells()). In \n"
    "other words, if reordering of the vertices does occur, the ordering of \n"
    "vertices in the array of corners will no longer refer to the same      \n"
    "triangulation.                                                         \n";



  const char generate_fixed_subdivided_parallelepiped_docstring[] =
    "A subdivided parallelepiped. The first corner point is the origin. The \n"
    "dim adjacent points are vectors describing the edges of the            \n"
    "parallelepiped with respect to the origin. Additional points are sums  \n"
    "of these dim vectors. The variable n_subdivisions designates the number\n"
    "of subdivisions in each of the dim directions. Colorizing is odne      \n"
    "according to hyper_rectangle().                                        \n";



  const char generate_varying_subdivided_parallelepiped_docstring[] =
    "A subdivided parallelepided, i.e., the same as above, but where the    \n"
    "number of subdivisions in each of the dim directions may vary.         \n"
    "Colorizing is done according to hyper_rectangle().                     \n";



  const char generate_enclosed_hyper_cube_docstring[] =
    "Hypercube with a layer of hypercubes around it. The first two          \n"
    "parameters give the lower and upper bound of the inner hypercube in all\n"
    "coordinate directions. thickness marks the size of the layer cells. If \n"
    "the flag colorize is set, the outer cells get material id's according  \n"
    "to the following scheme: extending over the inner cube (+/-)           \n"
    "x-direction: 1/2. In y-direction 4/8, in z-direction 16/32. The cells i\n"
    "at corners and edges (3d) get these values bitwise or'd.               \n";



  const char generate_hyper_ball_docstring[] =
    "Generate a hyperball, i.e. a circle or a ball around center with       \n"
    "given radius. In order to avoid degenerate cells at the boundaries,    \n"
    "the circle is triangulated by five cells, the ball by seven cells. The \n"
    "diameter of the center cell is chosen so that the aspect ratio of the  \n"
    "boundary cells after one refinement is optimized. You should attach a  \n"
    "SphericalManifold to the cells and faces for correct placement of      \n"
    "vertices upon refinement and to be able to use higher order mappings.  \n";



  const char generate_hyper_ball_balanced_docstring[] =
    "This is an alternative to hyper_ball with 12 cells in 2d and 32 cells  \n"
    "in 3d, which provides a better balance between the size of the cells   \n"
    "around the outer curved boundaries and the cell in the interior.       \n";



  const char generate_hyper_sphere_docstring[] =
    "Generate a hyper sphere, i.e., a surface of a ball in spacedim         \n"
    "dimensions. This function only exists for dim+1=spacedim in 2 and 3    \n"
    "space dimensions. You should attach a SphericalManifold to the cells   \n"
    "and faces for correct placement of vertices upon refinement and to be  \n"
    "able to use higher order mappings.                                     \n";



  const char generate_quarter_hyper_ball_docstring[] =
    "Generate a hyper-ball intersected with the positive orthant relate to  \n"
    "center, which contains three elements in 2d and four in 3d. The        \n"
    "boundary indicators for the final triangulations are 0 for the curved  \n"
    "boundary and 1 for the cut plane. The appropriate Manifold class is    \n"
    "SphericalManifold.                                                     \n";



  const char generate_hyper_shell_docstring[] =
    "Produce a hyper-shell, the region between two spheres around center,   \n"
    "with given inner_radius and outer_radius. The number n_cells indicates \n"
    "the number of cells of the resulting triangulation, i.e., how many     \n"
    "cells form the ring (in 2d) or the shell (in 3d).                      \n"
    "The appropriate manifold class is SphericalManifold.                   \n";



  const char generate_half_hyper_ball_docstring[] =
    "Generate a half hyper-ball around center, which contains four          \n"
    "elements in 2d and 6 in 3d. The cut plane is perpendicular to the      \n"
    "x-axis. The boundary indicators for the final triangulation are 0 for  \n"
    "the curved boundary and 1 for the cut plane. The appropriate manifold  \n"
    "class is SphericalManifold.                                            \n";



  const char generate_cylinder_docstring[] =
    "Create a dim dimensional cylinder where the x-axis serves as the axis \n"
    "of the cylinder. For the purposes of this function, a cylinder is     \n"
    "defined as a (dim - 1) dimensional disk of given radius, extruded     \n"
    "along the axis of the cylinder (which is the first coordinate         \n"
    "direction). Consequently, in three dimensions, the cylinder extends   \n"
    "from x=-half_length to x=+half_length and its projection into the     \n"
    "yz-plane is a circle of radius radius. In two dimensions, the         \n"
    "cylinder is a rectangle from x=-half_length to x=+half_length and     \n"
    "from y=-radius to y=radius.                                           \n";



  const char generate_subdivided_cylinder_docstring[] =
    "Create a dim dimensional cylinder where the x-axis serves as the axis \n"
    "of the cylinder. For the purposes of this function, a cylinder is     \n"
    "defined as a (dim - 1) dimensional disk of given radius, extruded     \n"
    "along the axis of the cylinder (which is the first coordinate         \n"
    "direction). Consequently, in three dimensions, the cylinder extends   \n"
    "from x=-half_length to x=+half_length and its projection into the     \n"
    "yz-plane is a circle of radius radius. In two dimensions, the         \n"
    "cylinder is a rectangle from x=-half_length to x=+half_length and     \n"
    "from y=-radius to y=radius. This function is only implemented for     \n"
    "dim == 3.                                                             \n";



  const char generate_truncated_cone_docstring[] =
    "Create a cut cone around the x-axis. The cone extends from            \n"
    "x=-half_length to x=half_length and its projection into the yz-plane  \n"
    "is a circle of radius radius_0 at x=-half_length and a circle of      \n"
    "radius radius_1 at x=+half_length. In between the radius is linearly  \n"
    "decreasing. In two dimensions, the cone is a trapezoid from           \n"
    "x=-half_length to x=+half_length and from y=-radius_0 to y=radius_0   \n"
    "at x=-half_length and from y=-radius_1 to y=radius_1 at               \n"
    "x=+half_length. In between the range of y is linearly decreasing.     \n";



  const char generate_hyper_cube_with_cylindrical_hole_docstring[] =
    "This function produces a square in the xy-plane with a cylindrical    \n"
    "hole in the middle. In 3d, this geometry is extruded in z direction   \n"
    "to the interval [0,L].                                                \n";



  const char generate_plate_with_a_hole_docstring[] =
    "Generate a rectangular plate with an (offset) cylindrical hole.       \n"
    "The geometry consists of 2 regions: The first is a square region      \n"
    "with length outer_radius and a hole of radius inner_radius.           \n"
    "The second region describes the remainder of the bulk material.       \n";



  const char generate_channel_with_cylinder_docstring[] =
    "Generate a grid consisting of a channel with a cylinder.              \n"
    "The channel has three distinct regions:                               \n"
    "  1. If n_shells is greater than zero, then there are that many       \n"
    "     shells centered around the cylinder,                             \n"
    "  2. a blending region between the shells and the rest of the         \n"
    "     triangulation, and                                               \n"
    "  3. a bulk region consisting of Cartesian cells.                     \n";



  const char shift_docstring[] =
    "Shift every vertex of the Triangulation by the given shift vector.     \n";



  const char scale_docstring[] =
    "Scale triangulation by a given scaling factor.                         \n";



  const char merge_triangulations_docstring[] =
    "Given two or more triangulations, create the triangulation that        \n"
    "contains the cells of given triangulations.                            \n";



  const char extrude_docstring[] =
    "Take a 2d Triangulation that is being extruded in z direction by       \n"
    "the total height of height using n_slices slices (minimum is 2).       \n"
    "The boundary indicators of the faces of input are going to be          \n"
    "assigned to the corresponding side walls in z direction. The           \n"
    "bottom and top get the next two free boundary indicators.              \n";



  const char flatten_triangulation_docstring[] =
    "Create a new flat triangulation out_tria which contains a single       \n"
    "level with all active cells of the input triangulation. If the spacedim\n"
    "are different, only the smalled spacedim components of the vertices are\n"
    "copied over. This is useful to create a Triangulation<2,3> out of a    \n"
    "Triangulation<2,2>, or to project a Triangulation<2,3> into a          \n"
    "Triangulation<2,2>, by neglecting the z component of the vertices. No  \n"
    "internal checks are performed on the vertices, which are assumed to    \n"
    "make sense topologically in the target spacedim dimensional space. If  \n"
    "this is not the case, you will encounter problems when using the       \n"
    "triangulation later on. All information about cell manifold_ids and    \n"
    "material ids are copied from one triangulation to the other, and only  \n"
    "the boundary manifold_ids and boundary_ids are copied over from the    \n"
    "faces of the triangulation to the faces of out_tria. If you need to    \n"
    "specify manifold ids on interior faces, they have to be specified      \n"
    "manually after the triangulation is created. This function will fail   \n"
    "if the input Triangulation contains hanging nodes.                     \n";



  const char replicate_docstring[] =
    "Replicate a given triangulation in multiple coordinate axes.          \n"
    "This function creates a new Triangulation equal to a dim-dimensional  \n"
    "array of copies of input.                                             \n";



  const char distort_random_docstring[] =
    "Distort the given triangulation by randomly moving around all the      \n"
    "vertices of the grid. The direction of movement of each vertex is      \n"
    "random, while the length of the shift vector has a value of factor     \n"
    "times the minimal length of the active edges adjacent to this vertex.  \n"
    "Note that factor should obviously be well below 0.5.                   \n";



  const char refine_global_docstring[] =
    "Refine all the cells times time.                                       \n";



  const char execute_coarsening_and_refinement_docstring[] =
    "Execute both refinement and coarsening of the Triangulation.           \n";



  const char active_cells_docstring[] =
    "Return the list of active cell accessors of the Triangulation.         \n";



  const char cells_docstring[] =
    "Return the list of cell accessors of the Triangulation.                \n";



  const char write_docstring[] =
    "Write the mesh to the output file according to the given data format.  \n"
    "The possible formats are:                                              \n"
    "  - none                                                               \n"
    "  - dx                                                                 \n"
    "  - gnuplot                                                            \n"
    "  - eps                                                                \n"
    "  - ucd                                                                \n"
    "  - xfig                                                               \n"
    "  - msh                                                                \n"
    "  - svg                                                                \n"
    "  - mathgl                                                             \n"
    "  - vtk                                                                \n"
    "  - vtu                                                                \n";



  const char read_docstring[] =
    "Read a mesh from the file according to the given data format.          \n"
    "The possible formats are:                                              \n"
    "  - msh                                                                \n"
    "  - vtk                                                                \n";



  const char save_docstring[] =
    "Write the Triangulation to a file.                                     \n";



  const char load_docstring[] =
    "Load the Triangulation from a file.                                    \n";



  const char set_manifold_docstring[] =
    "Assign a manifold object to a certain part of the triangulation.       \n"
    "The manifold_object is not copied and MUST persist until the           \n"
    "triangulation is destroyed.                                            \n";



  const char reset_manifold_docstring[] =
    "Reset those parts of the triangulation with the given manifold_number  \n"
    "to use a FlatManifold object.                                          \n";



  const char get_mesh_smoothing_docstring[] =
    "Return the mesh smoothing requirements that are obeyed.                \n";



  const char set_mesh_smoothing_docstring[] =
    "Set the mesh smoothing to mesh_smoothing.                              \n";



  const char transform_docstring[] =
    "Transform the vertices of the given triangulation by applying the      \n"
    "function object provided as first argument to all its vertices.        \n";



  const char find_active_cell_around_point_docstring[] =
    "Find and return an active cell that surrounds a given point p.         \n";



  const char compute_aspect_ratio_of_cells_docstring[] =
    "Computes an aspect ratio measure for all active cells and fills        \n"
    "a vector with one entry per cell.                                      \n";



  const char find_cells_adjacent_to_vertex_docstring[] =
    "Find and return a list of active cells that surround a given           \n"
    "vertex with index vertex_index.                                        \n";



  const char minimal_cell_diameter_docstring[] =
    "Return the diameter of the smallest active cell of a triangulation.    \n";



  const char maximal_cell_diameter_docstring[] =
    "Return the diameter of the largest active cell of a triangulation.    \n";



  const char convert_hypercube_to_simplex_mesh_docstring[] =
    "Convert a triangulation consisting only of hypercube cells            \n"
    "(quadrilaterals, hexahedra) to a triangulation only consisting of     \n"
    "simplices (triangles, tetrahedra).                                    \n";


  void
  export_triangulation()
  {
    boost::python::class_<TriangulationWrapper>(
      "Triangulation",
      boost::python::init<const std::string &,
                          boost::python::optional<const int, const bool>>(
        boost::python::args("dim",
                            "mesh_smoothing",
                            "check_for_distorted_cells")))
      .def(boost::python::init<const std::string &,
                               const std::string &,
                               boost::python::optional<const int, const bool>>(
        boost::python::args(
          "dim", "spacedim", "mesh_smoothing", "check_for_distorted_cells")))
      .def("n_active_cells",
           &TriangulationWrapper::n_active_cells,
           n_active_cells_docstring,
           boost::python::args("self"))
      .def("dim",
           &TriangulationWrapper::get_dim,
           dim_docstring,
           boost::python::args("self"))
      .def("spacedim",
           &TriangulationWrapper::get_spacedim,
           spacedim_docstring,
           boost::python::args("self"))
      .def("n_active_cells",
           &TriangulationWrapper::n_active_cells,
           n_active_cells_docstring,
           boost::python::args("self"))
      .def("n_cells",
           &TriangulationWrapper::n_cells,
           n_cells_docstring,
           boost::python::args("self"))
      .def("minimal_cell_diameter",
           &TriangulationWrapper::minimal_cell_diameter,
           minimal_cell_diameter_docstring,
           boost::python::args("self"))
      .def("maximal_cell_diameter",
           &TriangulationWrapper::maximal_cell_diameter,
           maximal_cell_diameter_docstring,
           boost::python::args("self"))
      .def("create_triangulation",
           &TriangulationWrapper::create_triangulation,
           create_triangulation_docstring,
           boost::python::args("self", "vertices", "cells_vertices"))
      .def("generate_hyper_cube",
           &TriangulationWrapper::generate_hyper_cube,
           generate_hyper_cube_overloads(
             boost::python::args("self", "left", "right", "colorize"),
             generate_hyper_cube_docstring))
      .def("generate_simplex",
           &TriangulationWrapper::generate_simplex,
           generate_simplex_docstring,
           boost::python::args("self", "vertices"))
      .def("generate_subdivided_hyper_cube",
           &TriangulationWrapper::generate_subdivided_hyper_cube,
           generate_subdivided_hyper_cube_overloads(
             boost::python::args("self", "repetitions", "left", "right"),
             generate_subdivided_hyper_cube_docstring))
      .def("generate_hyper_rectangle",
           &TriangulationWrapper::generate_hyper_rectangle,
           generate_hyper_rectangle_overloads(
             boost::python::args("self", "p1", "p2", "colorize"),
             generate_hyper_rectangle_docstring))
      .def("generate_subdivided_hyper_rectangle",
           &TriangulationWrapper::generate_subdivided_hyper_rectangle,
           generate_subdivided_hyper_rectangle_overloads(
             boost::python::args("self", "repetitions", "p1", "p2", "colorize"),
             generate_subdivided_hyper_rectangle_docstring))
      .def("generate_subdivided_steps_hyper_rectangle",
           &TriangulationWrapper::generate_subdivided_steps_hyper_rectangle,
           generate_subdivided_steps_hyper_rectangle_overloads(
             boost::python::args("self", "step_sizes", "p1", "p2", "colorize"),
             generate_subdivided_steps_hyper_rectangle_docstring))
      .def("generate_subdivided_material_hyper_rectangle",
           &TriangulationWrapper::generate_subdivided_material_hyper_rectangle,
           generate_subdivided_material_hyper_rectangle_overloads(
             boost::python::args(
               "self", "spacing", "p", "material_id", "colorize"),
             generate_subdivided_material_hyper_rectangle_docstring))
      .def("generate_hyper_cube_with_cylindrical_hole",
           &TriangulationWrapper::generate_hyper_cube_with_cylindrical_hole,
           generate_hyper_cube_with_cylindrical_hole_overloads(
             boost::python::args("self",
                                 "inner_radius",
                                 "outer_radius",
                                 "L",
                                 "repetitions",
                                 "colorize"),
             generate_hyper_cube_with_cylindrical_hole_docstring))
      .def("generate_cheese",
           &TriangulationWrapper::generate_cheese,
           generate_cheese_docstring,
           boost::python::args("self", "holes"))
      .def("generate_plate_with_a_hole",
           &TriangulationWrapper::generate_plate_with_a_hole,
           generate_plate_with_a_hole_overloads(
             boost::python::args("self",
                                 "inner_radius",
                                 "outer_radius",
                                 "pad_bottom",
                                 "pad_top",
                                 "pad_left",
                                 "pad_right",
                                 "center",
                                 "polar_manifold_id",
                                 "tfi_manifold_id",
                                 "L",
                                 "n_slices",
                                 "colorize"),
             generate_plate_with_a_hole_docstring))
      .def("generate_channel_with_cylinder",
           &TriangulationWrapper::generate_channel_with_cylinder,
           generate_channel_with_cylinder_overloads(
             boost::python::args(
               "shell_region_width", "n_shells", "skewness", "colorize"),
             generate_channel_with_cylinder_docstring))
      .def("generate_general_cell",
           &TriangulationWrapper::generate_general_cell,
           generate_general_cell_overloads(
             boost::python::args("self", "vertices", "colorize"),
             generate_general_cell_docstring))
      .def("generate_parallelogram",
           &TriangulationWrapper::generate_parallelogram,
           generate_parallelogram_overloads(
             boost::python::args("self", "corners", "colorize"),
             generate_parallelogram_docstring))
      .def("generate_parallelepiped",
           &TriangulationWrapper::generate_parallelepiped,
           generate_parallelepiped_overloads(
             boost::python::args("self", "corners", "colorize"),
             generate_parallelepiped_docstring))
      .def(
        "generate_fixed_subdivided_parallelepiped",
        &TriangulationWrapper::generate_fixed_subdivided_parallelepiped,
        generate_fixed_subdivided_parallelepiped_overloads(
          boost::python::args("self", "n_subdivisions", "corners", "colorize"),
          generate_fixed_subdivided_parallelepiped_docstring))
      .def(
        "generate_varying_subdivided_parallelepiped",
        &TriangulationWrapper::generate_varying_subdivided_parallelepiped,
        generate_varying_subdivided_parallelepiped_overloads(
          boost::python::args("self", "n_subdivisions", "corners", "colorize"),
          generate_varying_subdivided_parallelepiped_docstring))
      .def(
        "generate_enclosed_hyper_cube",
        &TriangulationWrapper::generate_enclosed_hyper_cube,
        generate_enclosed_hyper_cube_overloads(
          boost::python::args("self", "left", "right", "thickness", "colorize"),
          generate_enclosed_hyper_cube_docstring))
      .def("generate_hyper_ball",
           &TriangulationWrapper::generate_hyper_ball,
           generate_hyper_ball_overloads(
             boost::python::args("self", "center", "radius"),
             generate_hyper_ball_docstring))
      .def("generate_hyper_ball_balanced",
           &TriangulationWrapper::generate_hyper_ball_balanced,
           generate_hyper_ball_balanced_overloads(
             boost::python::args("self", "center", "radius"),
             generate_hyper_ball_balanced_docstring))
      .def("generate_hyper_sphere",
           &TriangulationWrapper::generate_hyper_sphere,
           generate_hyper_sphere_overloads(
             boost::python::args("self", "center", "radius"),
             generate_hyper_sphere_docstring))
      .def("generate_quarter_hyper_ball",
           &TriangulationWrapper::generate_quarter_hyper_ball,
           generate_quarter_hyper_ball_overloads(
             boost::python::args("self", "center", "radius"),
             generate_quarter_hyper_ball_docstring))
      .def("generate_half_hyper_ball",
           &TriangulationWrapper::generate_half_hyper_ball,
           generate_half_hyper_ball_overloads(
             boost::python::args("self", "center", "radius"),
             generate_half_hyper_ball_docstring))
      .def("generate_cylinder",
           &TriangulationWrapper::generate_cylinder,
           generate_cylinder_overloads(
             boost::python::args("self", "radius", "half_length"),
             generate_cylinder_docstring))
      .def("generate_subdivided_cylinder",
           &TriangulationWrapper::generate_subdivided_cylinder,
           generate_subdivided_cylinder_overloads(
             boost::python::args(
               "self", "x_subdivisions", "radius", "half_length"),
             generate_subdivided_cylinder_docstring))
      .def("generate_truncated_cone",
           &TriangulationWrapper::generate_truncated_cone,
           generate_truncated_cone_overloads(
             boost::python::args("self", "radius_0", "radius_2", "half_length"),
             generate_truncated_cone_docstring))
      .def("generate_hyper_shell",
           &TriangulationWrapper::generate_hyper_shell,
           generate_hyper_shell_overloads(boost::python::args("self",
                                                              "center",
                                                              "inner_radius",
                                                              "outer_radius",
                                                              "n_cells",
                                                              "colorize"),
                                          generate_hyper_shell_docstring))
      .def("shift",
           &TriangulationWrapper::shift,
           shift_docstring,
           boost::python::args("self", "shift"))
      .def("scale",
           &TriangulationWrapper::scale,
           scale_docstring,
           boost::python::args("self", "scaling_factor"))
      .def("merge_triangulations",
           &TriangulationWrapper::merge_triangulations,
           merge_triangulations_overloads(
             boost::python::args(
               "self", "triangulations", "vertex_tolerance", "copy_manifolds"),
             merge_triangulations_docstring))
      .def("extrude_triangulation",
           &TriangulationWrapper::extrude_triangulation,
           extrude_docstring,
           boost::python::args("self", "n_slices", "depth", "tria_out"))
      .def("flatten_triangulation",
           &TriangulationWrapper::flatten_triangulation,
           flatten_triangulation_docstring,
           boost::python::args("self", "tria_out"))
      .def("replicate_triangulation",
           &TriangulationWrapper::replicate_triangulation,
           replicate_docstring,
           boost::python::args("self", "tria_in", "extents"))
      .def("distort_random",
           &TriangulationWrapper::distort_random,
           distort_random_overloads(
             boost::python::args("self", "factor", "keep_boundary"),
             distort_random_docstring))
      .def("transform",
           &TriangulationWrapper::transform,
           transform_docstring,
           boost::python::args("self", "transformation"))
      .def("convert_hypercube_to_simplex_mesh",
           &TriangulationWrapper::convert_hypercube_to_simplex_mesh,
           convert_hypercube_to_simplex_mesh_docstring,
           boost::python::args("self", "tria_out"))
      .def("find_active_cell_around_point",
           &TriangulationWrapper::find_active_cell_around_point,
           find_active_cell_around_point_overloads(
             boost::python::args("self", "point", "mapping"),
             find_active_cell_around_point_docstring))
      .def("find_cells_adjacent_to_vertex",
           &TriangulationWrapper::find_cells_adjacent_to_vertex,
           find_cells_adjacent_to_vertex_docstring,
           boost::python::args("self", "vertex_index"))
      .def("compute_aspect_ratio_of_cells",
           &TriangulationWrapper::compute_aspect_ratio_of_cells,
           compute_aspect_ratio_of_cells_docstring,
           boost::python::args("self", "mapping", "quadrature"))
      .def("refine_global",
           &TriangulationWrapper::refine_global,
           refine_global_docstring,
           boost::python::args("self", "times"))
      .def("execute_coarsening_and_refinement",
           &TriangulationWrapper::execute_coarsening_and_refinement,
           execute_coarsening_and_refinement_docstring,
           boost::python::args("self"))
      .def("active_cells",
           &TriangulationWrapper::active_cells,
           active_cells_docstring,
           boost::python::args("self"))
      .def("cells",
           &TriangulationWrapper::cells,
           cells_docstring,
           boost::python::args("self"))
      .def("write",
           &TriangulationWrapper::write,
           write_docstring,
           boost::python::args("self", "filename", "format"))
      .def("read",
           &TriangulationWrapper::read,
           read_docstring,
           boost::python::args("self", "filename", "format"))
      .def("save",
           &TriangulationWrapper::save,
           save_docstring,
           boost::python::args("self", "filename"))
      .def("load",
           &TriangulationWrapper::load,
           load_docstring,
           boost::python::args("self", "filename"))
      .def("set_manifold",
           &TriangulationWrapper::set_manifold,
           set_manifold_docstring,
           boost::python::args("self", "number", "manifold"))
      .def("reset_manifold",
           &TriangulationWrapper::reset_manifold,
           reset_manifold_docstring,
           boost::python::args("self", "number"))
      .def("get_mesh_smoothing",
           &TriangulationWrapper::get_mesh_smoothing,
           get_mesh_smoothing_docstring,
           boost::python::args("self"))
      .def("set_mesh_smoothing",
           &TriangulationWrapper::set_mesh_smoothing,
           set_mesh_smoothing_docstring,
           boost::python::args("self", "mesh_smoothing"));

    boost::python::enum_<TriangulationWrapper::MeshSmoothing>("MeshSmoothing")
      .value("none", TriangulationWrapper::none)
      .value("limit_level_difference_at_vertices",
             TriangulationWrapper::limit_level_difference_at_vertices)
      .value("eliminate_unrefined_islands",
             TriangulationWrapper::eliminate_unrefined_islands)
      .value("patch_level_1", TriangulationWrapper::patch_level_1)
      .value("coarsest_level_1", TriangulationWrapper::coarsest_level_1)
      .value("allow_anisotropic_smoothing",
             TriangulationWrapper::allow_anisotropic_smoothing)
      .value("eliminate_refined_inner_islands",
             TriangulationWrapper::eliminate_refined_inner_islands)
      .value("eliminate_refined_boundary_islands",
             TriangulationWrapper::eliminate_refined_boundary_islands)
      .value("do_not_produce_unrefined_islands",
             TriangulationWrapper::do_not_produce_unrefined_islands)
      .value("smoothing_on_refinement",
             TriangulationWrapper::smoothing_on_refinement)
      .value("smoothing_on_coarsening",
             TriangulationWrapper::smoothing_on_coarsening)
      .value("maximum_smoothing", TriangulationWrapper::maximum_smoothing);
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
