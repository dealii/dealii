// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test distributed solution_transfer with one processor

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "coarse_grid_common.h"

template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(){};

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double f = p[0] * 2.0 + 1.0;
    if (dim > 1)
      f *= p[1] * 3.3 - 1.0;
    if (dim > 2)
      f *= p[2] * 5.0;
    return f;
  };
};

template <int dim>
void
test(std::ostream & /*out*/)
{
  MyFunction<dim>                           func;
  MappingQ<dim>                             mapping(1);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(3);
  DoFHandler<dim>        dofh(tr);
  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs(fe);

  SolutionTransfer<dim, Vector<double>> soltrans(dofh);

  for (int i = 0; i < 4; ++i)
    {
      for (typename Triangulation<dim>::cell_iterator cell = tr.begin();
           cell != tr.end();
           ++cell)
        {
          if (cell->has_children() && !(Testing::rand() % 3))
            for (unsigned int c = 0;
                 c < GeometryInfo<dim>::max_children_per_cell;
                 ++c)
              if (!cell->child(c)->has_children())
                cell->child(c)->set_coarsen_flag();
        }
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        {
          if (!(Testing::rand() % 3))
            {
              cell->clear_coarsen_flag();
              cell->set_refine_flag();
            }
        }

      tr.prepare_coarsening_and_refinement();

      Vector<double> solution(dofh.n_dofs());
      VectorTools::interpolate(mapping, dofh, func, solution);

      soltrans.prepare_for_coarsening_and_refinement(solution);

      tr.execute_coarsening_and_refinement();
      if (0)
        {
          std::ofstream out("grid-1.eps");
          GridOut       grid_out;
          grid_out.write_eps(tr, out);
        }

      dofh.distribute_dofs(fe);

      Vector<double> interpolated_solution(dofh.n_dofs());
      soltrans.interpolate(interpolated_solution);

      deallog << "norm: " << interpolated_solution.l2_norm() << std::endl;
      Vector<double> difference(tr.n_global_active_cells());

      VectorTools::integrate_difference(mapping,
                                        dofh,
                                        interpolated_solution,
                                        func,
                                        difference,
                                        QGauss<dim>(2),
                                        VectorTools::L2_norm);
      deallog << "error: " << difference.l2_norm() << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  initlog();
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#else
  (void)argc;
  (void)argv;
#endif

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
