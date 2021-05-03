// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// common framework to check whether an element of polynomial order p can
// represent functions of order q for projection from quadrature points.

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>

#include "../tests.h"

template <int dim>
class F : public Function<dim>
{
public:
  F(const unsigned int q, const unsigned int n_components)
    : Function<dim>(n_components)
    , q(q)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    Assert((component == 0) && (this->n_components == 1), ExcInternalError());
    double val = 0;
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int i = 0; i <= q; ++i)
        val += (d + 1) * (i + 1) * std::pow(p[d], 1. * i);
    return val;
  }

  VectorizedArray<double>
  value(const Point<dim, VectorizedArray<double>> &p_vec) const
  {
    VectorizedArray<double> res = make_vectorized_array(0.);
    Point<dim>              p;
    for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
      {
        for (unsigned int d = 0; d < dim; d++)
          p[d] = p_vec[d][v];
        res[v] = value(p);
      }
    return res;
  }


private:
  const unsigned int q;
};

struct QData
{
  double density;
};


template <typename VectorType, int dim>
void
do_project(const parallel::distributed::Triangulation<dim> &triangulation,
           const FiniteElement<dim> &                       fe,
           const unsigned int                               p)
{
  AssertThrow(fe.n_components() == 1, ExcNotImplemented());
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  deallog << "n_cells=" << triangulation.n_global_active_cells() << std::endl;
  deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  for (unsigned int q = 0; q <= p; ++q)
    {
      // setup quadrature data:
      F<dim>     function(q, fe.n_components());
      VectorType field_relevant;
      double     field_l2_norm = 0.;
      {
        QGauss<dim> quadrature_formula(p + 2);
        CellDataStorage<typename parallel::distributed::
                          Triangulation<dim, dim>::cell_iterator,
                        QData>
          qp_manager;
        {
          FEValues<dim> fe_values(fe,
                                  quadrature_formula,
                                  update_quadrature_points);

          const unsigned int n_q_points = quadrature_formula.size();

          typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                qp_manager.initialize(cell, quadrature_formula.size());
                std::vector<std::shared_ptr<QData>> qpd =
                  qp_manager.get_data(cell);
                for (unsigned int q = 0; q < n_q_points; ++q)
                  qpd[q]->density =
                    function.value(fe_values.quadrature_point(q));
              }
        }

        VectorType field(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
        VectorTools::project<dim, VectorType, dim>(
          MappingQGeneric<dim>(1),
          dof_handler,
          constraints,
          quadrature_formula,
          [=](const typename DoFHandler<dim>::active_cell_iterator &cell,
              const unsigned int q) -> double {
            return qp_manager.get_data(cell)[q]->density;
          },
          field);

        field_relevant.reinit(dof_handler.locally_owned_dofs(),
                              locally_relevant_dofs,
                              MPI_COMM_WORLD);
        field_relevant = field;

        field_l2_norm = field.l2_norm();
      }

      // L2 norm of the difference between FE field and the function
      double L2_norm = 0.;
      {
        QGauss<dim>   quadrature_formula_error(std::max(p, q) + 1);
        FEValues<dim> fe_values(fe,
                                quadrature_formula_error,
                                update_values | update_quadrature_points |
                                  update_JxW_values);

        const unsigned int  dofs_per_cell = fe.dofs_per_cell;
        const unsigned int  n_q_points    = quadrature_formula_error.size();
        std::vector<double> values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler
                                                                .begin_active(),
                                                       endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);
              fe_values.get_function_values(field_relevant, values);

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                L2_norm +=
                  Utilities::fixed_power<2>(
                    values[q_point] -
                    function.value(fe_values.quadrature_point(q_point))) *
                  fe_values.JxW(q_point);
            }
      }

      L2_norm = Utilities::MPI::sum(L2_norm, MPI_COMM_WORLD);
      L2_norm = std::sqrt(L2_norm);

      deallog << fe.get_name() << ", P_" << q
              << ", rel. error=" << L2_norm / field_l2_norm << std::endl;

      if (q <= p)
        if (L2_norm > 1e-10 * field_l2_norm)
          deallog << "Projection failed with relative error "
                  << L2_norm / field_l2_norm << std::endl;
    }
}



// check the given element of polynomial order p. the last parameter, if
// given, denotes a gap in convergence order; for example, the Nedelec element
// of polynomial degree p has normal components of degree p-1 and therefore
// can only represent polynomials of degree p-1 exactly. the gap is then 1.
template <typename VectorType, int dim>
void
test_no_hanging_nodes(const FiniteElement<dim> &fe, const unsigned int p)
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  do_project<VectorType>(triangulation, fe, p);
}



// same test as above, but this time with a mesh that has hanging nodes
template <typename VectorType, int dim>
void
test_with_hanging_nodes(const FiniteElement<dim> &fe, const unsigned int p)
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  if (triangulation.begin_active()->is_locally_owned())
    triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global(1);

  do_project<VectorType>(triangulation, fe, p);
}
