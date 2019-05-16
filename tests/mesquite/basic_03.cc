/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * A demonstration of how to hand-construct a mesh to be smoothed.
 * Note: The mesh coordinates are smoothed "in-place" and are not copied!
 *
 * This simple program is part of a tutorial found in the user guide.
 *
 * Knupp, P.; Freitag-Diachin, L. & Tidwell, B.
 * Mesquite Mesh Quality Improvement Toolkit User's Guide
 * Sandia National Laboratories, Sandia National Laboratories, 2013
 *
 * See sec 4.2 Accessing Mesh In Arrays
 */


#include <deal.II/grid/grid_tools_mesquite.h>

#include <Mesquite_all_headers.hpp>

#include <fstream>

#include "../tests.h"


void
test_2d()
{
  // Error tracker
  Mesquite2::MsqError err;

  // ---------------------

  // Mesh vertex coordinates
  // Vertex coordinates must be stored in an array of double-precision
  // floating-point values. The coor- dinate values must be interleaved,
  // meaning that the x, y, and z coordinate values for a single vertex
  // are contiguous in memory.
  double coords[] = {0, 0, 1, 0, 2, 0, 0, 1, 0.5, 0.5, 2, 1, 0, 2, 1, 2, 2, 2};

  // Quadrilateral element connectivity
  // The element connectivity (vertices in each element) must be stored
  // in an interleaved format as an array of long integers. The vertices
  // in each element are specified by an integer i, where the location of
  // the coordinates of the corresponding vertex is located at position
  // 3*i in the vertex coordinates array.
  const unsigned long quads[] = {
    0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7};

  // Constraints: All vertices except for the middle one are fixed.
  // The fixed / boundary state of the vertices must be stored in an array
  // of integer values, where a value of 1 indicates a fixed vertex and
  // a value of 0 indicates a free vertex. The values in this array must
  // be in the same order as the corresponding vertex coordinates in the
  // coordinate array.
  const int fixed[] = {1, 1, 1, 1, 0, 1, 1, 1, 1};

  // Pass the above mesh into Mesquite
  // The mesh must be composed of a single element type.
  // NOTE: When using the ArrayMesh interface, the application is responsible
  // for
  //       managing the storage of the mesh data. The ArrayMesh does NOT copy
  //       the input mesh.
  Mesquite2::ArrayMesh mesh(
    2,      // specify a 2-d mesh (two coordinate values per vertex)
    9,      // number of vertices
    coords, // vertex coordinates
    fixed,  // constraint flags
    4,      // number of elements
    Mesquite2::QUADRILATERAL, // element type
    quads                     // element connectivity
  );

  // ---------------------

  // Surface to constrain the 2d elements to
  Mesquite2::PlanarDomain domain(Mesquite2::PlanarDomain::XY);

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &domain);

  // Output the old location of the center vertex
  deallog << "Old vertex location: ( " << coords[8] << ", " << coords[9] << " )"
          << std::endl;

  // Improve the mesh
  Mesquite2::ShapeImprover shape_wrapper;
  shape_wrapper.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Output the new location of the center vertex
  deallog << "New vertex location: ( " << coords[8] << ", " << coords[9] << " )"
          << std::endl;
}


void
test_3d()
{
  // Error tracker
  Mesquite2::MsqError err;

  // ---------------------

  // Mesh vertex coordinates
  // Vertex coordinates must be stored in an array of double-precision
  // floating-point values. The coordinate values must be interleaved,
  // meaning that the x, y, and z coordinate values for a single vertex
  // are contiguous in memory.
  double coords[] = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0.5, 0.5,
                     0, 2, 1, 0, 0, 2, 0, 1, 2, 0, 2, 2, 0};

  // Quadrilateral element connectivity
  // The element connectivity (vertices in each element) must be stored
  // in an interleaved format as an array of long integers. The vertices
  // in each element are specified by an integer i, where the location of
  // the coordinates of the corresponding vertex is located at position
  // 3*i in the vertex coordinates array.
  const unsigned long quads[] = {
    0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7};

  // Constraints: All vertices except for the middle one are fixed.
  // The fixed / boundary state of the vertices must be stored in an array
  // of integer values, where a value of 1 indicates a fixed vertex and
  // a value of 0 indicates a free vertex. The values in this array must
  // be in the same order as the corresponding vertex coordinates in the
  // coordinate array.
  const int fixed[] = {1, 1, 1, 1, 0, 1, 1, 1, 1};

  // Pass the above mesh into Mesquite
  // The mesh must be composed of a single element type.
  // NOTE: When using the ArrayMesh interface, the application is responsible
  // for
  //       managing the storage of the mesh data. The ArrayMesh does NOT copy
  //       the input mesh.
  Mesquite2::ArrayMesh mesh(
    3,      // specify a 3-d mesh (three coordinate values per vertex)
    9,      // number of vertices
    coords, // vertex coordinates
    fixed,  // constraint flags
    4,      // number of elements
    Mesquite2::QUADRILATERAL, // element type
    quads                     // element connectivity
  );

  // ---------------------

  // Surface to constrain the 2d elements to
  Mesquite2::PlanarDomain domain(Mesquite2::PlanarDomain::XY);

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &domain);

  // Output the old location of the center vertex
  deallog << "Old vertex location: ( " << coords[12] << ", " << coords[13]
          << ", " << coords[14] << " )" << std::endl;

  // Improve the mesh
  Mesquite2::ShapeImprover shape_wrapper;
  shape_wrapper.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Output the new location of the center vertex
  deallog << "New vertex location: ( " << coords[12] << ", " << coords[13]
          << ", " << coords[14] << " )" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  {
    test_2d();
  }
  deallog.pop();

  deallog.push("3d");
  {
    test_3d();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
