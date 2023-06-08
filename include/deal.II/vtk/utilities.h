// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef dealii_vtk_utilities_h
#define dealii_vtk_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

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

#endif
#endif
