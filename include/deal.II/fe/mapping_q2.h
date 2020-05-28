// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_mapping_q2_h
#define dealii_mapping_q2_h

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Implementation of a quadratic mapping from the reference cell to a general
 * quadrilateral/hexahedron given 9-noded quadrangles or 27-noded hexahedra.
 *
 * The mapping implemented by this class is similar to creating a mapping of
 * the type MappingQ(2) or MappingQGeneric(2). However, the difference is that
 * this mapping gets the necessary support points from the input file read
 * with read_msh() instead of computing them or getting them from a manifold.
 *
 * @note This class does not support mesh refinement yet. Also, only Gmsh
 * input files can currently be read and processed to give the necessary
 * vector of support points.
 *
 * @author Daniel Paukner, Peter Munch 2020
 */

template <int dim, int spacedim = dim>
class MappingQ2 : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] support_points
   * A vector of a vector of points that stores the support points associated
   * with a cell. Note that only the support points are stored in this vector,
   * not the vertices.
   */
  MappingQ2(std::vector<std::vector<Point<spacedim>>> &support_points);

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;


private:
  /**
   * Returns the support points necessary for the quadratic mapping. However,
   * instead of computing the support points, this function just returns the
   * support points stored in the support_points vector corresponding to the
   * current cell.
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * Vector with the support points read from the mesh file.
   */
  std::vector<std::vector<Point<spacedim>>> support_points;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
