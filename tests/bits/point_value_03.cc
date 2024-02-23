// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check VectorTools::point_value
// This test is the parallel version of point_value_01



#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


using VectorType = LinearAlgebra::distributed::Vector<double>;

template <int dim>
class MySquareFunction : public Function<dim, VectorType::value_type>
{
public:
  MySquareFunction()
    : Function<dim, VectorType::value_type>()
  {}

  virtual VectorType::value_type
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return (component + 1) * p.square() + 1;
  }

  virtual void
  vector_value(const Point<dim>               &p,
               Vector<VectorType::value_type> &values) const
  {
    values(0) = value(p, 0);
  }
};


template <int dim>
class MyExpFunction : public Function<dim>
{
public:
  MyExpFunction()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return std::exp(p(0));
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    values(0) = value(p, 0);
  }
};



template <int dim>
void
make_mesh(Triangulation<dim> &tria)
{
  GridGenerator::hyper_cube(tria, -1, 1);

  // refine the mesh in a random way so as to
  // generate as many cells with
  // hanging nodes as possible
  tria.refine_global(4 - dim);
  const double steps[4] = {/*d=0*/ 0, 7, 3, 3};
  for (unsigned int i = 0; i < steps[dim]; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
        tria.begin_active();
      for (unsigned int index = 0; cell != tria.end(); ++cell, ++index)
        if (index % (3 * dim) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
}



template <int dim>
void
check()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  make_mesh(tria);

  FE_Q<dim>       element(QIterated<1>(QTrapezoid<1>(), 3));
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(element);

  static const MySquareFunction<dim> function;

  VectorType v;

  const IndexSet locally_owned_dofs = dof.locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof);

  v.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  VectorTools::interpolate(dof, function, v);
  v.update_ghost_values();

  // v = 1.;

  // for the following points, check the
  // function value, output it, and
  // verify that the value retrieved from
  // the interpolated function is close
  // enough to that of the real function
  //
  // also verify that the actual value is
  // roughly correct
  Point<dim> p[3];
  for (unsigned int d = 0; d < dim; ++d)
    {
      p[0][d] = 0;
      p[1][d] = 0.5;
      p[2][d] = 1. / 3.;
    }
  Vector<VectorType::value_type> value(1);
  const auto                    &mapping = get_default_linear_mapping(tria);
  for (unsigned int i = 0; i < 3; ++i)
    {
      {
        bool point_in_locally_owned_cell = false;
        value(0)                         = 0;
        {
          auto cell_and_ref_point =
            GridTools::find_active_cell_around_point(mapping, dof, p[i]);
          if (cell_and_ref_point.first.state() == IteratorState::valid)
            {
              point_in_locally_owned_cell =
                cell_and_ref_point.first->is_locally_owned();
            }
        }
        if (point_in_locally_owned_cell)
          {
            VectorTools::point_value(mapping, dof, v, p[i], value);
          }
      }
      value(0) = Utilities::MPI::sum(value(0), MPI_COMM_WORLD);
      deallog << -value(0) << std::endl;
      std::cout << value(0) << std::endl;

      Assert(std::abs(value(0) - function.value(p[i])) < 2e-4,
             ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;


  MPILogInitAll log;

  deallog << std::setprecision(4);

  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
