// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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


/**
 * @defgroup simplex Simplex support (experimental)
 *
 * @brief This module describes the experimental simplex support in deal.II.
 *
 * Simplex and mixed meshes in deal.II are still experimental, i.e., work
 * in progress. Large parts of the library have been ported to be able to
 * operate on such kind of meshes. However, there are still many functions
 * that need to be generalized. You can get a good overview of the ported
 * functionalities by taking a look at the tests in the folder
 * "tests/simplex". In the following, we provide two very basic examples
 * to get started and provide some implementation details.
 *
 * @section simplex_reference_example_simplex Example: simplex mesh
 *
 * The following code shows how to work with simplex meshes:
 *
 * @include step_3_simplex.cc
 *
 * @section simplex_reference_example_mixed Example: mixed mesh
 *
 * The following code shows how to work with mixed meshes:
 *
 * @include step_3_mixed.cc
 *
 * @section simplex_reference_cells Reference cells
 *
 * In 2D, we provide triangles and quadrilaterals with the following possible
 * orientations in 3D:
 *
 * <div class="twocolumn" style="width: 100%">
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_0.png
 *     </div>
 *     <div class="text" align="center">
 *       2D: triangle and quadrilateral
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html reference_cells_1.png
 *     </div>
 *     <div class="text" align="center">
 *       Possible orientations of triangles and quadrilaterals in 3D
 *     </div>
 *   </div>
 * </div>
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
