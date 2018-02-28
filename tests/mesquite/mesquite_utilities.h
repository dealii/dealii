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
 * A collection of utilities for the Mesquite tests
 */

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <Mesquite_all_headers.hpp>

#include <fstream>

#include "../tests.h"


std::string
translate_error_code(const Mesquite2::MsqError::ErrorCode &code)
{
  // The meanings to the error codes can be found in the file
  // <path_to_mesquite_source>/src/Misc/Mesquite_MsqError.hpp
  // or at
  // https://github.com/trilinos/mesquite/blob/master/src/Misc/Mesquite_MsqError.hpp#L109
  switch (code)
    {
      case (Mesquite2::MsqError::NO_ERROR):
        return "NO_ERROR: No error";
        break;
      case (Mesquite2::MsqError::UNKNOWN_ERROR):
        return "UNKNOWN_ERROR: Unknown error occured";
        break;
      case (Mesquite2::MsqError::OUT_OF_MEMORY):
        return "OUT_OF_MEMORY: Unable to allocate the necessary memory";
        break;
      case (Mesquite2::MsqError::INVALID_ARG):
        return "INVALID_ARG: Invalid function argument passed";
        break;
      case (Mesquite2::MsqError::NOT_INITIALIZED):
        return "NOT_INITIALIZED: Object not initialized";
        break;
      case (Mesquite2::MsqError::INVALID_STATE):
        return "INVALID_STATE: Object is in an invalid state";
        break;
      case (Mesquite2::MsqError::FILE_ACCESS):
        return "FILE_ACCESS: File cannot be opened/created";
        break;
      case (Mesquite2::MsqError::FILE_FORMAT):
        return "FILE_FORMAT: Wrong file type";
        break;
      case (Mesquite2::MsqError::PARSE_ERROR):
        return "PARSE_ERROR: Error parsing input (or input file)";
        break;
      case (Mesquite2::MsqError::IO_ERROR):
        return "IO_ERROR: An I/O error occured (e.g. read from file failed.)";
        break;
      case (Mesquite2::MsqError::INVALID_MESH):
        return "INVALID_MESH: The mesh is invalid";
        break;
      case (Mesquite2::MsqError::NO_PD_STORAGE_MODE):
        return "NO_PD_STORAGE_MODE: No storage mode chosen within PatchData";
        break;
      case (Mesquite2::MsqError::NOT_IMPLEMENTED):
        return "NOT_IMPLEMENTED: Requested functionality is not (yet) implemented";
        break;
      case (Mesquite2::MsqError::INTERNAL_ERROR):
        return "INTERNAL_ERROR: A bug in Mesquite";
        break;
      case (Mesquite2::MsqError::INTERRUPTED):
        return "INTERRUPTED: Application or user interrupted operation";
        break;
      case (Mesquite2::MsqError::TAG_ALREADY_EXISTS):
        return "TAG_ALREADY_EXISTS: Attempt to create tag that already exists";
        break;
      case (Mesquite2::MsqError::TAG_NOT_FOUND):
        return "TAG_NOT_FOUND: Specified tag does not exist";
        break;
      case (Mesquite2::MsqError::UNSUPPORTED_ELEMENT):
        return "UNSUPPORTED_ELEMENT: The element type is not supported.";
        break;
      case (Mesquite2::MsqError::PARALLEL_ERROR):
        return "PARALLEL_ERROR: An error occurred in parallel";
        break;
      case (Mesquite2::MsqError::BARRIER_VIOLATED):
        return "BARRIER_VIOLATED: Barrier violated when processing barrier Target Metric";
        break;
      default:
        AssertThrow(code < Mesquite2::MsqError::LAST_ERROR_CODE,
                    ExcMessage("Invalid error code"));
        break;
    }
  return "";
}

std::string
translate_error_code(const Mesquite2::MsqError &err)
{
  return translate_error_code(err.error_code());
}


/**
 * An error message for Mesquite functions.
 */
DeclException1(ExcMesquiteError,
               Mesquite2::MsqError,
               << "An error with error number " << arg1.error_code()
               << " occurred while calling a Mesquite function."
               << "\n    " << translate_error_code(arg1));


void smooth_mesquite_in_place_arraymesh(Triangulation<2> &tria)
{
  const unsigned int dim = 2;

  // ---------------------
  // Convert triangulation to a format that Mesquite understands

  // Mesh vertex coordinates
  // Note: This has a stride of dim between vertices
  std::vector<double>            coords(dim * tria.n_vertices());
  const std::vector<Point<dim>> &vertices = tria.get_vertices();
  for (unsigned int v = 0; v < vertices.size(); ++v)
    for (unsigned int d = 0; d < dim; ++d)
      {
        const unsigned int c = dim * v + d; // Flattened coordinate index
        Assert(c < coords.size(), ExcIndexRange(c, 0, coords.size()));
        coords[c] = vertices[v][d];
      }

  // Element connectivity
  // Note: This has a stride of n_vertices_per_cell between elements
  std::vector<unsigned long> cells(GeometryInfo<dim>::vertices_per_cell *
                                   tria.n_active_cells());
  unsigned int               cell_count = 0;
  for (auto cell : tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          const unsigned int c =
            GeometryInfo<dim>::vertices_per_cell * cell_count +
            v; // Flattened global vertex index
          Assert(c < cells.size(), ExcIndexRange(c, 0, cells.size()));
          cells[c] = cell->vertex_index(v);
        }
      ++cell_count;
    }

  // Constraints: All boundary vertices are fixed.
  // Note: This has a stride of 1 between vertices
  std::vector<int> constraints(tria.n_vertices(), 0);
  for (auto cell : tria.active_cell_iterators())
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
          if (cell->face(f)->at_boundary() == false)
            continue;
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
               ++v)
            {
              const unsigned int c = cell->face(f)->vertex_index(v);
              Assert(c < constraints.size(),
                     ExcIndexRange(c, 0, constraints.size()));
              constraints[c] = 1;
            }
        }
    }

  // Pass the above into Mesquite
  // Note: This Mesquite2::ArrayMesh object understands both
  //       2d and 3d entities, while a Mesquite2::MeshImpl can
  //       only work in 3d.
  Mesquite2::ArrayMesh mesh(
    dim,                   // specify a mesh (dim coordinate values per vertex)
    tria.n_vertices(),     // number of vertices
    coords.data(),         // vertex coordinates
    constraints.data(),    // constraint flags
    tria.n_active_cells(), // number of elements
    Mesquite2::QUADRILATERAL, // element type
    cells.data()              // element connectivity
  );

  // ---------------------

  // Error tracker
  Mesquite2::MsqError err;

  // Surface to constrain the 2d elements to
  Mesquite2::PlanarDomain domain(Mesquite2::PlanarDomain::XY);

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &domain);

  // Improve the mesh
  Mesquite2::ShapeImprover shape_wrapper;
  shape_wrapper.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // ---------------------
  // Move vertices in the triangulation as required
  std::vector<bool> vertex_touched(tria.n_vertices(), false);
  cell_count = 0;
  for (auto cell : tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          const unsigned int idx = cell->vertex_index(v);
          Assert(idx < vertex_touched.size(),
                 ExcIndexRange(idx, 0, vertex_touched.size()));
          if (!vertex_touched[idx])
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const unsigned int c =
                    dim * idx + d; // Flattened coordinate index
                  Assert(c < coords.size(), ExcIndexRange(c, 0, coords.size()));
                  cell->vertex(v)[d] = coords[c];
                }

              vertex_touched[idx] = true;
            }
        }
      ++cell_count;
    }
}


template <int dim>
void
output_mesh(const Triangulation<dim> &triangulation,
            std::string               filename,
            const bool                write_to_deallog = true)
{
  GridOut grid_out;

  if (dim == 2)
    {
      filename += ".eps";
      std::ofstream out(filename);
      grid_out.write_eps(triangulation, out);

      if (write_to_deallog)
        grid_out.write_eps(triangulation, deallog.get_file_stream());
    }
  else
    {
      filename += ".vtk";
      std::ofstream out(filename);
      grid_out.write_vtk(triangulation, out);

      if (write_to_deallog)
        grid_out.write_vtk(triangulation, deallog.get_file_stream());
    }
}
