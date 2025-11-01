// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// document bug in the new AffineConstraints<double>::distribute():
/*
--------------------------------------------------------
An error occurred in line <1028> of file
</scratch/deal-trunk/deal.II/include/deal.II/lac/constraint_matrix.templates.h>
in function void dealii::AffineConstraints<double>::distribute(VectorType&)
const [with VectorType = dealii::TrilinosWrappers::MPI::Vector] The violated
condition was: vec(it->line) == it->entries.size() The name and call sequence of
the exception was: ExcIncorrectConstraint(it->line, it->entries.size())
Additional Information:
While distributing the constraint for DoF 41, it turns out that one of the
processors who own the 2 degrees of freedom that x_41 is constrained against
does not know about the constraint on x_41. Did you not initialize the
AffineConstraints<double> with the appropriate locally_relevant set so that
every processor who owns a DoF that constrains another DoF also knows about this
constraint?
--------------------------------------------------------

 */

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  if (myid == 0)
    {
      typename Triangulation<dim>::active_cell_iterator it = tr.begin_active();
      ++it;
      ++it;
      ++it;
      for (unsigned int i = 0; i < 5; ++i)
        {
          it->set_refine_flag();
          ++it;
        }
    }
  tr.execute_coarsening_and_refinement();

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  const IndexSet locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dofh);

  AffineConstraints<double> cm;
  cm.reinit(dofh.locally_owned_dofs(), locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, cm);

  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  Functions::ZeroFunction<dim> homogeneous_dirichlet_bc(1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  VectorTools::interpolate_boundary_values(dofh, dirichlet_boundary, cm);

  cm.close();

  TrilinosWrappers::MPI::Vector vec(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  deallog << "locally owned:" << std::endl;
  vec.locally_owned_elements().print(deallog);
  deallog << "relevant set:" << std::endl;
  locally_relevant_set.print(deallog);
  if (myid == 0)
    {
      deallog << "constraint_matrix:" << std::endl;
      cm.print(deallog.get_file_stream());
    }

  cm.distribute(vec);

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test<2>();
    }
  else
    test<2>();
}
