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
                             FiniteElementData<dim>::H1),
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



template <int dim>
void
FE_BernardiRaugel<dim>::fill_fe_values(
  const typename Triangulation<dim, dim>::cell_iterator &cell,
  const CellSimilarity::Similarity                       cell_similarity,
  const Quadrature<dim> &                                quadrature,
  const Mapping<dim, dim> &                              mapping,
  const typename Mapping<dim, dim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, dim>
    &                                                       mapping_data,
  const typename FiniteElement<dim, dim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
    &output_data) const
{
  std::cout << "fill_fe_values called!" << std::endl;
  std::cout << "Calling FE_PolyTensor<dim, dim>::fill_fe_values" << std::endl;

  // Call the base implementation then we'll touch up the data slightly
  FE_PolyTensor<dim, dim>::fill_fe_values(
     cell,
     cell_similarity,
     quadrature,
     mapping,
     mapping_internal,
     mapping_data,
     fe_internal,
     output_data);
  std::cout << "Length of mapping_data.normal_vectors" << normal_vectors.size() << std::endl;

  // Convert to the correct internal data class for this FE class.
  // Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
  //    ExcInternalError());
  const typename FE_PolyTensor<dim>::InternalData &fe_data =
    static_cast<const typename FE_PolyTensor<dim>::InternalData &>(fe_internal);

  const unsigned int n_q_points = quadrature.size();

  // convert the output data to use the outward normal from the physical cell
  // i.e. if the original shape function is \hat{\psi} = \vec{\hat{n}}*phi,
  // then we recover psi and compute
  // \psi = \vec{n}*(\hat{\psi}_1^2 + \hat{\psi}_2^2 + \hat{\psi}_3^2)^(1/2)
  std::cout << "FE_BR adjusting bubbles:" << std::endl;
  for (unsigned int i = dim*GeometryInfo<dim>::vertices_per_cell; i < this->dofs_per_cell; ++i)
  {
    const unsigned int first =
            output_data.shape_function_to_row_table[i * this->n_components() +
                                                    this->get_nonzero_components(i)
                                                      .first_selected_component()];

    for (unsigned int k = 0; k < n_q_points; ++k)
    {
      double psi = 0;
      unsigned int f = i - dim*GeometryInfo<dim>::vertices_per_cell;
      Tensor<1,dim> normal = (*cell)->face(f)->normal_vector(k);

      for (unsigned int d = 0; d < dim; ++d)
        psi += fe_data.shape_values[i][k][d]*fe_data.shape_values[i][k][d];

      psi = sqrt(psi);

      for (unsigned int d = 0; d < dim; ++d)
      {
        output_data.shape_values(first + d, k) =
          normal[d]*psi;
      }

    }
  }

  std::cout << "Done!" << std::endl;
}

template class FE_BernardiRaugel<2>;
template class FE_BernardiRaugel<3>;

DEAL_II_NAMESPACE_CLOSE
