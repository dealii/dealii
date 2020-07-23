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


// Define an FE_Q1_Nonlocal finite element, where dofs are not distributed
// on vertices (with dpo[0]=1) but on cells (with dpo[dim+2]=vertices_per_cell)
// and the non-local dof numeration is based on the vertex id

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

#include <iostream>

#include "../tests.h"

template <int dim>
std::vector<unsigned int>
get_dpo_vector()
{
  std::vector<unsigned int> dpo(dim + 2, 0.0);
  dpo[dim + 1] = GeometryInfo<dim>::vertices_per_cell;
  return dpo;
}


/**
 * Associates non local dofs according to the vertex index.
 */
template <int dim>
class NonLocalQ1DoFHandler : public NonLocalDoFHandler<dim>
{
public:
  virtual std::vector<types::global_dof_index>
  get_non_local_dof_indices(
    const DoFCellAccessor<dim, dim, false> &accessor) const override
  {
    if (!tria)
      tria = &(accessor.get_triangulation());
    std::vector<types::global_dof_index> dofs(accessor.n_vertices());
    for (unsigned int i = 0; i < dofs.size(); ++i)
      dofs[i] = accessor.vertex_index(i);
    return dofs;
  }

  virtual types::global_dof_index
  n_additional_non_local_dofs() const override
  {
    Assert(tria, ExcInternalError());
    return tria->n_vertices();
  }

private:
  mutable SmartPointer<const Triangulation<dim>> tria = nullptr;
};



template <int dim>
class FE_Q1_Nonlocal : public FE_Q_Base<TensorProductPolynomials<dim>, dim, dim>
{
public:
  FE_Q1_Nonlocal(const std::shared_ptr<NonLocalQ1DoFHandler<dim>> &ptr =
                   std::shared_ptr<NonLocalQ1DoFHandler<dim>>())
    : FE_Q_Base<TensorProductPolynomials<dim>, dim, dim>(
        TensorProductPolynomials<dim>(
          Polynomials::generate_complete_Lagrange_basis(
            QGaussLobatto<1>(2).get_points())),
        FiniteElementData<dim>(get_dpo_vector<dim>(),
                               1,
                               1,
                               FiniteElementData<dim>::H1),
        std::vector<bool>(1, false))
    , non_local_dh(ptr ? ptr : std::make_shared<NonLocalQ1DoFHandler<dim>>())
  {
    this->unit_support_points = QTrapez<dim>().get_points();
  }

  virtual std::unique_ptr<FiniteElement<dim>>
  clone() const override
  {
    return std::make_unique<FE_Q1_Nonlocal<dim>>(non_local_dh);
  }

  virtual std::string
  get_name() const override
  {
    return "FE_Q_Nonlocal<dim>";
  }

  virtual std::shared_ptr<const NonLocalDoFHandler<dim>>
  get_non_local_dof_handler() const override
  {
    return non_local_dh;
  }

  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override
  {
    AssertDimension(support_point_values.size(),
                    this->get_unit_support_points().size());
    AssertDimension(support_point_values.size(), nodal_values.size());
    AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      {
        AssertDimension(support_point_values[i].size(), 1);
        nodal_values[i] = support_point_values[i](0);
      }
  }

private:
  std::shared_ptr<NonLocalQ1DoFHandler<dim>> non_local_dh;
};
