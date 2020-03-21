// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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
// represent functions of order q for projection from quadrature points of
// matrix-free approach.

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

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
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>

#include "../tests.h"

/**
 * A function of polynomial degree @p q with n_components.
 */
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


/**
 * A general function to perform projection of polynomial function of degrees
 * [0,p]
 * onto the FE space given by @p fe and @p triangulation .
 */
template <int fe_degree, int n_q_points_1d, int dim>
void
do_project(const parallel::distributed::Triangulation<dim> &triangulation,
           const std::vector<const FiniteElement<dim> *> &  fes,
           const unsigned int                               p,
           const unsigned int                               fe_index = 0)
{
  // we only project scalar value function
  AssertThrow(fes[fe_index]->n_components() == 1, ExcNotImplemented());

  deallog << "n_cells=" << triangulation.n_global_active_cells() << std::endl;

  std::vector<std::shared_ptr<DoFHandler<dim>>> dof_handlers(fes.size());
  std::vector<IndexSet>                  locally_relevant_dofs(fes.size());
  std::vector<AffineConstraints<double>> constraints(fes.size());

  std::vector<const DoFHandler<dim> *>           dof_handlers_mf(fes.size());
  std::vector<const AffineConstraints<double> *> constraints_mf(fes.size());
  for (unsigned int i = 0; i < fes.size(); ++i)
    {
      dof_handlers[i] = std::make_shared<DoFHandler<dim>>(triangulation);
      dof_handlers[i]->distribute_dofs(*fes[i]);
      deallog << "n_dofs=" << dof_handlers[i]->n_dofs() << std::endl;

      DoFTools::extract_locally_relevant_dofs(*dof_handlers[i],
                                              locally_relevant_dofs[i]);

      constraints[i].reinit(locally_relevant_dofs[i]);
      DoFTools::make_hanging_node_constraints(*dof_handlers[i], constraints[i]);
      constraints[i].close();

      dof_handlers_mf[i] = dof_handlers[i].get();
      constraints_mf[i]  = &constraints[i];
    }

  QGauss<1> quadrature_formula_1d(n_q_points_1d);

  Table<2, VectorizedArray<double>> qp_data;

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::partition_color;
  additional_data.mapping_update_flags =
    update_values | update_JxW_values | update_quadrature_points;
  std::shared_ptr<MatrixFree<dim, double>> data(new MatrixFree<dim, double>());
  data->reinit(dof_handlers_mf,
               constraints_mf,
               quadrature_formula_1d,
               additional_data);

  for (unsigned int q = 0; q <= p; ++q)
    {
      // setup quadrature data:
      F<dim> function(q, fes[fe_index]->n_components());

      // initialize a quadrature data
      {
        FEEvaluation<dim, fe_degree, n_q_points_1d, 1, double> fe_eval(
          *data, fe_index);
        const unsigned int n_cells    = data->n_macro_cells();
        const unsigned int n_q_points = fe_eval.n_q_points;

        qp_data.reinit(n_cells, n_q_points);
        for (unsigned int cell = 0; cell < n_cells; ++cell)
          {
            fe_eval.reinit(cell);
            for (unsigned int q = 0; q < n_q_points; ++q)
              qp_data(cell, q) = function.value(fe_eval.quadrature_point(q));
          }
      }

      LinearAlgebra::distributed::Vector<double> field;
      data->initialize_dof_vector(field, fe_index);
      VectorTools::project<dim, LinearAlgebra::distributed::Vector<double>>(
        data,
        constraints[fe_index],
        n_q_points_1d,
        [=](const unsigned int cell, const unsigned int q)
          -> VectorizedArray<double> { return qp_data(cell, q); },
        field,
        fe_index);

      field.update_ghost_values();

      const double field_l2_norm = field.l2_norm();

      // L2 norm of the difference between FE field and the function
      double L2_norm = 0.;
      {
        QGauss<dim>   quadrature_formula_error(std::max(p, q) + 1);
        FEValues<dim> fe_values(*fes[fe_index],
                                quadrature_formula_error,
                                update_values | update_quadrature_points |
                                  update_JxW_values);

        const unsigned int  dofs_per_cell = fes[fe_index]->dofs_per_cell;
        const unsigned int  n_q_points    = quadrature_formula_error.size();
        std::vector<double> values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handlers[fe_index]->begin_active(),
          endc = dof_handlers[fe_index]->end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);
              fe_values.get_function_values(field, values);

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

      deallog << fes[fe_index]->get_name() << ", P_" << q
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
template <int fe_degree, int n_q_points_1d, int dim>
void
test_no_hanging_nodes(const FiniteElement<dim> &fe, const unsigned int p)
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  do_project<fe_degree, n_q_points_1d, dim>(triangulation, {{&fe}}, p);
}



// same test as above, but this time with a mesh that has hanging nodes
template <int fe_degree, int n_q_points_1d, int dim>
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

  do_project<fe_degree, n_q_points_1d, dim>(triangulation, {{&fe}}, p);
}


// same as above but for multiple fes
template <int fe_degree, int n_q_points_1d, int dim>
void
test_with_hanging_nodes(const std::vector<const FiniteElement<dim> *> &fes,
                        const unsigned int                             p,
                        const unsigned int                             fe_index)
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  if (triangulation.begin_active()->is_locally_owned())
    triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global(1);

  do_project<fe_degree, n_q_points_1d, dim>(triangulation, fes, p, fe_index);
}
