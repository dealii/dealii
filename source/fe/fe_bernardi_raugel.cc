// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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


#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernardi_raugel.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>
#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

template <int dim>
FE_BernardiRaugel<dim>::FE_BernardiRaugel(const unsigned int p)
  : FE_PolyTensor<dim>(
      PolynomialsBernardiRaugel<dim>(p),
      FiniteElementData<dim>(get_dpo_vector(),
                             dim,
                             2,
                             FiniteElementData<dim>::Hdiv),
      std::vector<bool>(PolynomialsBernardiRaugel<dim>::n_polynomials(p), true),
      std::vector<ComponentMask>(PolynomialsBernardiRaugel<dim>::n_polynomials(
                                   p),
                                 std::vector<bool>(dim, true)))
{
  Assert(dim == 2 || dim == 3, ExcImpossibleInDim(dim));
  Assert(p == 1, ExcMessage("Only BR1 elements are available"));

  // const unsigned int n_dofs = this->dofs_per_cell;

  this->mapping_kind = {mapping_none};
  // These must be done first, since
  // they change the evaluation of
  // basis functions

  // Set up the generalized support
  // points
  initialize_support_points();
}



template <int dim>
std::string
FE_BernardiRaugel<dim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_BR<" << dim << ">(" << 1 << ")";

  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_BernardiRaugel<dim>::clone() const
{
  return std::make_unique<FE_BernardiRaugel<dim>>(*this);
}



template <int dim>
void
FE_BernardiRaugel<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  AssertDimension(support_point_values[0].size(), dim);
  Assert(nodal_values.size() == this->dofs_per_cell,
         ExcDimensionMismatch(nodal_values.size(), this->dofs_per_cell));

  std::vector<Tensor<1, dim>> normals;
  for (unsigned int i : GeometryInfo<dim>::face_indices())
    {
      Tensor<1, dim> normal;
      normal[i / 2] = 1;
      normals.push_back(normal);
    }

  for (unsigned int i = 0; i < dim * GeometryInfo<dim>::vertices_per_cell; ++i)
    nodal_values[i] = support_point_values[i][i % dim];

  // Compute the support points values for the bubble functions
  for (unsigned int i = dim * GeometryInfo<dim>::vertices_per_cell;
       i < dim * GeometryInfo<dim>::vertices_per_cell +
             GeometryInfo<dim>::faces_per_cell;
       ++i)
    {
      nodal_values[i] = 0;
      for (unsigned int j = 0; j < dim; ++j)
        nodal_values[i] +=
          support_point_values[i][j] *
          normals[i - dim * GeometryInfo<dim>::vertices_per_cell][j];
    }
}



template <int dim>
std::vector<unsigned int>
FE_BernardiRaugel<dim>::get_dpo_vector()
{
  // compute the number of unknowns per cell interior/face/edge
  //
  // there are <tt>dim</tt> degrees of freedom per vertex and there
  // is 1 degree of freedom per edge in 2D (face in 3D)
  std::vector<unsigned int> dpo(dim + 1, 0u);
  dpo[0]       = dim;
  dpo[dim - 1] = 1u;

  return dpo;
}



template <int dim>
void
FE_BernardiRaugel<dim>::initialize_support_points()
{
  // The support points for our shape functions are the vertices and
  // the face midpoints, for a total of #vertices + #faces points
  this->generalized_support_points.resize(this->dofs_per_cell);

  // We need dim copies of each vertex for the first dim*vertices_per_cell
  // generalized support points
  for (unsigned int i = 0; i < dim * GeometryInfo<dim>::vertices_per_cell; ++i)
    this->generalized_support_points[i] =
      GeometryInfo<dim>::unit_cell_vertex(i / dim);

  // The remaining 2*dim points are the edge midpoints
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < 2; ++j)
        {
          Point<dim> p;
          p[0] = 0.5;
          p[1] = 0.5;
          if (dim == 3)
            p[2] = 0.5;
          p[i] = j;

          const unsigned int k =
            dim * GeometryInfo<dim>::vertices_per_cell + i * 2 + j;
          this->generalized_support_points[k] = p;
        }
    }
}

template class FE_BernardiRaugel<2>;
template class FE_BernardiRaugel<3>;

DEAL_II_NAMESPACE_CLOSE
