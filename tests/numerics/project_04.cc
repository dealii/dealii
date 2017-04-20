// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check VectorTools::Project for quadrature data


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/quadrature_point_data.h>

#include <fstream>

namespace LA
{
  using namespace ::LinearAlgebraTrilinos;
}

struct QData
{
  double density;
};


// define the multi-linear function x or x*y or x*y*z that we will
// subsequently project onto the ansatz space
template <int dim>
class F : public Function<dim>
{
public:
  virtual double value (const Point<dim> &p,
                        const unsigned int = 0) const
  {
    double s = 1;
    for (unsigned int i=0; i<dim; ++i)
      s *= p[i];
    return s;
  }
};


template<int dim>
void test()
{
  F<dim> function;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss<dim> quadrature_formula(3);

  CellDataStorage<typename parallel::distributed::Triangulation<dim,dim>::cell_iterator,QData> qp_manager;
  {
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_quadrature_points);

    const unsigned int   n_q_points    = quadrature_formula.size();

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          qp_manager.initialize(cell, quadrature_formula.size());
          std::vector<std::shared_ptr<QData> > qpd = qp_manager.get_data(cell);
          for (unsigned int q=0; q<n_q_points; ++q)
            qpd[q]->density = function.value(fe_values.quadrature_point(q));
        }
  }

  ConstraintMatrix constraints;
  constraints.close ();

  VectorTools::Projector<LA::MPI::Vector,
              LA::MPI::SparseMatrix,
              LA::SolverCG,
              LA::MPI::PreconditionAMG> projector;
  projector.initialize(dof_handler,
                       constraints,
                       quadrature_formula,
                       MPI_COMM_WORLD);

  LA::MPI::Vector field(dof_handler.locally_owned_dofs(),
                        MPI_COMM_WORLD);
  projector.project<dim,dim>(dof_handler,
                             constraints,
                             quadrature_formula,
                             MPI_COMM_WORLD,
                             [=] (const typename DoFHandler<dim>::active_cell_iterator & cell, const unsigned int q) -> double { return qp_manager.get_data(cell)[q]->density; },
                             field);

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           locally_relevant_dofs);

  LA::MPI::Vector field_relevant(dof_handler.locally_owned_dofs(),
                                 locally_relevant_dofs,
                                 MPI_COMM_WORLD);
  field_relevant = field;



  // L2 norm of the difference between FE field and the function
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values  | update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  std::vector<double>  values (n_q_points);

  double L2_norm = 0.;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(field_relevant,
                                      values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          L2_norm += Utilities::fixed_power<2>(values[q_point] -  function.value(fe_values.quadrature_point(q_point))) *
                     fe_values.JxW(q_point);

      }

  L2_norm = Utilities::MPI::sum(L2_norm, MPI_COMM_WORLD);
  L2_norm = std::sqrt(L2_norm);
  deallog << L2_norm << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}

