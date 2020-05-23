// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_ball_overloads,
                                         generate_hyper_ball,
                                         1,
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

  const char n_active_cells_docstring[] =
    "Return the number of active cells.                                     \n";



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
    "lexicographic ordering (see GridReordering::reoder_grid()). In other   \n"
    "words, if reordering of the vertices does occur, the ordering of       \n"
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



  const char generate_hyper_cube_with_cylindrical_hole_docstring[] =
    "This function produces a square in the xy-plane with a cylindrical     \n"
    "hole in the middle. In 3d, this geometry is extruded in z direction    \n"
    "to the interval [0,L].                                                 \n";



  const char shift_docstring[] =
    "Shift every vertex of the Triangulation by the given shift vector.     \n";



  const char scale_docstring[] =
    "Scale triangulation by a given scaling factor.                         \n";



  const char merge_docstring[] =
    "Given two triangulations, create the triangulation that contains       \n"
    "the cells of both triangulations.                                      \n";



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


  void
  export_triangulation()
  {
    boost::python::class_<TriangulationWrapper>(
      "Triangulation",
      boost::python::init<const std::string &>(boost::python::args("dim")))
      .def(boost::python::init<const std::string &, const std::string &>(
        boost::python::args("dim", "spacedim")))
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
           merge_docstring,
           boost::python::args("self", "triangulation_1", "triangulation_2"))
      .def("extrude_triangulation",
           &TriangulationWrapper::extrude_triangulation,
           extrude_docstring,
           boost::python::args("self", "n_slices", "depth", "tria_out"))
      .def("flatten_triangulation",
           &TriangulationWrapper::flatten_triangulation,
           flatten_triangulation_docstring,
           boost::python::args("self", "tria_out"))
      .def("distort_random",
           &TriangulationWrapper::distort_random,
           distort_random_overloads(
             boost::python::args("self", "factor", "keep_boundary"),
             distort_random_docstring))
      .def("transform",
           &TriangulationWrapper::transform,
           transform_docstring,
           boost::python::args("self", "transformation"))
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
           boost::python::args("self", "number"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
