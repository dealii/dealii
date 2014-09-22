// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// Test that DofRenumbering::hierarchical() gives the same DoFs independent of
// the number of processors. We do this by constructing a vector where each
// element depends on the location and the index. We output the l2_norm() of
// that vector and hope the output is that same.
// This test fails right now, because p4est gets the coarse cells in a different
// order.


#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <boost/concept_check.hpp>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr (MPI_COMM_WORLD);

  GridGenerator::hyper_shell (tr,
                              Point<dim>(),
                              0.3,
                              1.0,
                              12, true);
//  GridGenerator::hyper_cube (tr, -1.0, 1.0);
  tr.refine_global (1);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      if (cell->center()[0] < 0.0)
        {
          cell->set_refine_flag();
        }

  tr.execute_coarsening_and_refinement ();

  DoFHandler<dim> dofh (tr);

  static const FE_Q<dim> fe (1);
  dofh.distribute_dofs (fe);
  DoFRenumbering::hierarchical (dofh);

  if (myid == 0)
    deallog << "n_global_active_cells: " << tr.n_global_active_cells() << std::endl;

  TrilinosWrappers::MPI::Vector vector;
  vector.reinit (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  {
    ConstraintMatrix cm;
    cm.close();
    QGauss<dim>   quadrature_formula (2);
    FEValues<dim> fe_values (fe, quadrature_formula, update_quadrature_points  |
                             update_JxW_values | update_values);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    Vector<double> local_vector (dofs_per_cell);
    const unsigned int   n_q_points      = quadrature_formula.size();

    typename DoFHandler<dim>::active_cell_iterator
    cell = dofh.begin_active(),
    endc = dofh.end();
    for (; cell != endc; ++cell)
      if (cell->subdomain_id() == tr.locally_owned_subdomain())
        {
          fe_values.reinit (cell);
          cell->get_dof_indices (local_dof_indices);
          local_vector = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  local_vector (i) += fe_values.shape_value (i, q) *
                                      1.0 * local_dof_indices[i] *
                                      (fe_values.quadrature_point (q) [0] +
                                       fe_values.quadrature_point (q).square())
                                      * fe_values.JxW (q);
                }
            }
          cm.distribute_local_to_global (local_vector, local_dof_indices, vector);
        }
    vector.compress(VectorOperation::add);
  }
  double norm = vector.l2_norm();
  if (myid == 0)
    deallog << "Norm: " << norm << std::endl;
}


int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push (Utilities::int_to_string (myid));

  if (myid == 0)
    {
      std::ofstream logfile ("output");
      deallog.attach (logfile);
      deallog.depth_console (0);
      deallog.threshold_double (1.e-10);

      deallog.push ("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
