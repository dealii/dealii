// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vtk_utilities_h
#define dealii_vtk_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>

#ifdef DEAL_II_WITH_VTK

#  include <vtkDoubleArray.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Interface to the Visualization Toolkit (VTK).
 *
 * VTK is an open-source, freely available software system for 3D computer
 * graphics, modeling, image processing, volume rendering, scientific
 * visualization, and 2D plotting.

 * It supports a wide variety of visualization algorithms and advanced
 * modeling techniques, and it takes advantage of both threaded and distributed
 * memory parallel processing for speed and scalability, respectively.
 *
 * You can learn more about the VTK library at https://vtk.org/
 */
namespace VTKWrappers
{
  /**
   * Convert from a deal.II Point to a VTK double array.
   *
   * @tparam dim Dimension of the point
   * @param [in] p An input deal.II Point<dim>
   * @return A VTK smart pointer to the data array.
   */
  template <int dim>
  inline vtkSmartPointer<vtkDoubleArray>
  dealii_point_to_vtk_array(const dealii::Point<dim> &p);

  /**
   * @brief Read a VTK mesh file and populate a deal.II Triangulation.
   *
   * This function reads the mesh from the specified VTK file and fills the
   * given Triangulation object. If cleanup is true, overlapping points in the
   * VTK file are merged using VTK's cleaning utilities.
   *
   * @param vtk_filename The name of the input VTK file.
   * @param tria The Triangulation object to populate.
   * @param cleanup If true, merge overlapping points in the VTK file (default: true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  template <int dim, int spacedim>
  void
  read_tria(const std::string            &vtk_filename,
            Triangulation<dim, spacedim> &tria,
            const bool                    cleanup            = true,
            const double                  relative_tolerance = 0.0);

  /**
   * Read cell data (scalar or vector) from a VTK file and store it in the
   * output vector.
   *
   * This function reads the specified cell data array (scalar or vector) from
   * the given VTK file and stores it in the provided output vector. For vector
   * data, all components are stored in row-major order (cell0_comp0,
   * cell0_comp1, ..., cell1_comp0, ...).
   *
   * @param vtk_filename The name of the input VTK file.
   * @param cell_data_name The name of the cell data array to read.
   * @param output_vector The vector to store the cell data values.
   */
  void
  read_cell_data(const std::string &vtk_filename,
                 const std::string &cell_data_name,
                 Vector<double>    &output_vector);

  /**
   * Read vertex data from a VTK file and store it in the output vector.
   *
   * This function reads the specified vertex data array (scalar or vector) from
   * the given VTK file and stores it in the provided output vector.
   *
   * For vector data, all components are stored in row-major order
   * (vertex0_comp0, vertex0_comp1, ..., vertex1_comp0, ...).
   *
   * @param vtk_filename The name of the input VTK file.
   * @param vertex_data_name The name of the vertex data array to read.
   * @param output_vector The vector to store the vertex data values.
   */
  void
  read_vertex_data(const std::string &vtk_filename,
                   const std::string &vertex_data_name,
                   Vector<double>    &output_vector);

  /**
   * Create a FiniteElement representation for data stored in a VTK file.
   *
   * This function inspects the specified VTK file and builds a FiniteElement
   * (typically an FESystem) that can store the data fields contained in the
   * file. It returns a pair consisting of a std::unique_ptr to the constructed
   * FiniteElement and a vector of field names corresponding to the VTK data
   * arrays.
   *
   * Only meshes with a single cell type are supported (either all quads/hexes
   * or all triangles/tetrahedra). Mixed meshes are not supported. The mesh type
   * is inferred from the first cell encountered in the VTK file.
   *
   * Mapping rules:
   * - Vertex data (VTK point data):
   *   - quad/hex meshes: represented by FE_Q(1) blocks;
   *   - simplex meshes: represented by FE_SimplexP(1) blocks.
   * - Cell data (VTK cell data):
   *   - quad/hex meshes: represented by FE_DGQ(0) blocks;
   *   - simplex meshes: represented by FE_SimplexDGP(0) blocks.
   *
   *   Vector-valued cell or vertex data are represented as an FESystem of the
   *   scalar cell element repeated for each component.
   *
   * @param vtk_filename The path to the input VTK file.
   */
  template <int dim, int spacedim>
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
            std::vector<std::string>>
  vtk_to_finite_element(const std::string &vtk_filename);

#  ifndef DOXYGEN
  // Template implementations

  template <int dim>
  inline vtkSmartPointer<vtkDoubleArray>
  dealii_point_to_vtk_array(const dealii::Point<dim> &p)
  {
    vtkSmartPointer<vtkDoubleArray> p_vtk =
      vtkSmartPointer<vtkDoubleArray>::New();

    p_vtk->SetNumberOfComponents(dim);
    p_vtk->SetNumberOfTuples(1);

    for (int d = 0; d < dim; ++d)
      p_vtk->FillComponent(d, p[d]);

    return p_vtk;
  }

#  endif

} // namespace VTKWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
