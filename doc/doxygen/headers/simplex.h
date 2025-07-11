// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @defgroup simplex Simplex support
 *
 * This group describes the simplex support in deal.II.
 *
 * Simplex and mixed meshes in deal.II are work in progress.
 * Large parts of the library have been ported to be able to
 * operate on such kind of meshes. However, there are still many functions
 * that need to be generalized.
 *
 * @section simplex_functionality_list Important Simplex Functionality
 *
 * Here is an incomplete list of functionality related to simplex
 * computations:
 * - Mesh generation:
 *   GridGenerator::implicit_function(),
 *   GridGenerator::convert_hypercube_to_simplex_mesh(),
 *   GridGenerator::subdivided_hyper_rectangle_with_simplices(),
 *   GridGenerator::subdivided_hyper_cube_with_simplices()
 * - Quadratures:
 *   QGaussWedge, QGaussSimplex, QWitherdenVincentSimplex
 * - FiniteElements:
 *   FE_SimplexP, FE_SimplexDGP, FE_SimplexP_Bubbles
 *   FE_PyramidP, FE_PyramidDGP, FE_WedgeP, FE_WedgeDGP
 * - Mapping:
 *   MappingFE
 * - Other:
 *   GridIn::read_vtk(), GridIn::read_msh(), GridIn::read_comsol_mphtxt()
 *
 *
 *
 * @section Examples
 *
 * You can get a good overview of the ported
 * functionalities by taking a look at the tests in the folder
 * "tests/simplex". In the following, we provide two very basic examples
 * to get you started and to provide some implementation details.
 *
 * @subsection simplex_reference_example_simplex Example: simplex mesh
 *
 * The following code shows how to work with simplex meshes:
 *
 * @include step_3_simplex.cc
 *
 * @subsection simplex_reference_example_mixed Example: mixed mesh
 *
 * The following code shows how to work with mixed meshes:
 *
 * @include step_3_mixed.cc
 *
 * @section simplex_reference_cells Reference cells
 *
 * In 2D, we provide triangles and quadrilaterals:
 *
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_0.png
 *     </div>
 *     <div class="text" align="center">
 *       2D: triangle and quadrilateral
 *     </div>
 *   </div>
 *
 * In 3D, tetrahedra, pyramids, wedges, and hexahedra are available:
 *
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_2.png
 *     </div>
 *     <div class="text" align="center">
 *       3D: Tetrahedron
 *     </div>
 *   </div>
 *
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_3.png
 *     </div>
 *     <div class="text" align="center">
 *       3D: Pyramid
 *     </div>
 *   </div>
 *
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_4.png
 *     </div>
 *     <div class="text" align="center">
 *       3D: Wedge
 *     </div>
 *   </div>
 *
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_5.png
 *     </div>
 *     <div class="text" align="center">
 *       3D: Hexahedron
 *     </div>
 *   </div>
 *
 * Each surface of a 3D reference cell consists of 2D reference cells. The
 * documentation of the enumeration of the numbering of their vertices and
 * lines are given in the right columns.
 *
 */
