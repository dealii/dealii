// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test regression in FEFieldFunction that doesn't throw
// VectorTools::ExcPointNotAvailableHere()

/*
An error occurred in line <3405> of file
</ssd/deal-git/include/deal.II/dofs/dof_accessor.templates.h> in function void
dealii::DoFCellAccessor<dealii::DoFHandler<2, 2>, false>::get_dof_values(const
InputVector &, ForwardIterator, ForwardIterator) const [DoFHandlerType =
dealii::DoFHandler<2, 2>, lda = false, InputVector =
dealii::TrilinosWrappers::MPI::Vector, ForwardIterator = double *] The violated
condition was: this->is_artificial() == false Additional information: Can't ask
for DoF indices on artificial cells.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class LinearFunction : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p[0] + 2;
  }
};


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube_with_cylindrical_hole(tr,
                                                  0.05 /* cylinder radius */,
                                                  0.4 / 2.0 /* box radius */);
  tr.refine_global(1);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
                                             MPI_COMM_WORLD);
  VectorTools::interpolate(dofh, LinearFunction<dim>(), interpolated);

  const IndexSet relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = interpolated;

  typename Functions::FEFieldFunction<dim, TrilinosWrappers::MPI::Vector>

    field_func(dofh, x_rel);

  Point<2>              p(0.1, 0.0);
  std::vector<Point<2>> points;
  points.push_back(p);

  for (int i = 0;; ++i)
    {
      try
        {
          if (i == 0)
            {
              deallog << "vector_value:" << std::endl;
              Vector<double> out(1);
              field_func.vector_value(p, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 1)
            {
              deallog << "value:" << std::endl;
              double out = field_func.value(p);
              deallog << "  OK: " << out << std::endl;
            }
          else if (i == 2)
            {
              deallog << "value_list:" << std::endl;
              std::vector<double> out(1, -42.0f);
              field_func.value_list(points, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 3)
            {
              deallog << "vector_value_list:" << std::endl;
              std::vector<Vector<double>> out(1, Vector<double>(1));
              field_func.vector_value_list(points, out);
              deallog << "  OK: " << out[0][0] << std::endl;
            }
          else if (i == 4)
            {
              deallog << "vector_gradient:" << std::endl;
              std::vector<Tensor<1, 2>> out(1, Tensor<1, 2>());
              field_func.vector_gradient(p, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 5)
            {
              deallog << "gradient:" << std::endl;
              Tensor<1, 2> out = field_func.gradient(p);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 6)
            {
              deallog << "vector_gradient_list:" << std::endl;
              std::vector<std::vector<Tensor<1, 2>>> out(
                1, std::vector<Tensor<1, 2>>(1, Tensor<1, 2>()));
              field_func.vector_gradient_list(points, out);
              deallog << "  OK: " << out[0][0] << std::endl;
            }
          else if (i == 7)
            {
              deallog << "gradient_list:" << std::endl;
              std::vector<Tensor<1, 2>> out(1, Tensor<1, 2>());
              field_func.gradient_list(points, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 8)
            {
              deallog << "laplacian:" << std::endl;
              double out = field_func.laplacian(p);
              deallog << "  OK: " << out << std::endl;
            }
          else if (i == 9)
            {
              deallog << "vector_laplacian:" << std::endl;
              Vector<double> out(1);
              field_func.vector_laplacian(p, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 10)
            {
              deallog << "laplacian_list:" << std::endl;
              std::vector<double> out(1);
              field_func.laplacian_list(points, out);
              deallog << "  OK: " << out[0] << std::endl;
            }
          else if (i == 11)
            {
              deallog << "vector_laplacian_list:" << std::endl;
              std::vector<Vector<double>> out(1, Vector<double>(1));
              field_func.vector_laplacian_list(points, out);
              deallog << "  OK: " << out[0][0] << std::endl;
            }
          else
            break;
        }
      catch (...)
        {
          deallog << "  ExcPointNotAvailableHere" << std::endl;
        }
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  deal_II_exceptions::disable_abort_on_exception();

  test<2>();
}
