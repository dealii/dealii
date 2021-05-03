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



// Check CellDataStorage class. To that end first evaluate some
// quadratic function at quadrature points. Then refine cells and project using
// FE_Q(2). Finally check that the values at quadrature points are still
// consistent with the original function.


#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(1)
  {}

  double
  value(const Point<dim> &p, const unsigned int comp = 0) const
  {
    const double x = p[0];
    const double y = p[1];
    // some function we know we can project with FE_Q<dim>(2)
    return 0.5 * x * x + 2.1 * y * y + 2;
  }
};

class MyQData : public TransferableQuadraturePointData
{
public:
  MyQData(){};
  virtual ~MyQData(){};

  unsigned int
  number_of_values() const
  {
    return 1;
  }

  double value;

  virtual void
  pack_values(std::vector<double> &scalars) const
  {
    Assert(scalars.size() == 1, ExcInternalError());
    scalars[0] = value;
  }

  virtual void
  unpack_values(const std::vector<double> &scalars)
  {
    Assert(scalars.size() == 1, ExcInternalError());
    value = scalars[0];
  }
};

const double eps = 1e-10;
DeclException3(ExcWrongValue,
               double,
               double,
               double,
               << arg1 << " != " << arg2 << " with delta = " << arg3);


/**
 * Loop over quadrature points and check that value is the same as given by the
 * function.
 */
template <int dim, typename DATA>
void
check_qph(parallel::distributed::Triangulation<dim> &tr,
          const CellDataStorage<typename Triangulation<dim, dim>::cell_iterator,
                                DATA> &              manager,
          const Quadrature<dim> &                    rhs_quadrature,
          const MyFunction<dim> &                    func)
{
  DoFHandler<dim> dof_handler(tr);
  FE_Q<dim>       dummy_fe(1);
  FEValues<dim>   fe_values(dummy_fe, rhs_quadrature, update_quadrature_points);
  dof_handler.distribute_dofs(dummy_fe);
  typename Triangulation<dim, dim>::active_cell_iterator cell;
  for (cell = tr.begin_active(); cell != tr.end(); ++cell)
    if (cell->is_locally_owned())
      {
        typename DoFHandler<dim>::active_cell_iterator dof_cell(*cell,
                                                                &dof_handler);
        fe_values.reinit(dof_cell);
        const std::vector<Point<dim>> &q_points =
          fe_values.get_quadrature_points();
        const std::vector<std::shared_ptr<const DATA>> qpd =
          manager.get_data(cell);
        for (unsigned int q = 0; q < q_points.size(); q++)
          {
            const double value  = func.value(q_points[q]);
            const double value2 = qpd[q]->value;
            AssertThrow(std::fabs(value - value2) < eps,
                        ExcWrongValue(value, value2, value - value2));
          }
      }
  dof_handler.clear();
}

template <int dim>
void
test()
{
  unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  const MyFunction<dim> my_func;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tr, 2);
  tr.refine_global(1);
  typename Triangulation<dim, dim>::active_cell_iterator cell;

  // pppulate quadrature point data
  QGauss<dim> rhs(4);
  CellDataStorage<typename Triangulation<dim, dim>::cell_iterator, MyQData>
    data_storage;
  parallel::distributed::ContinuousQuadratureDataTransfer<dim, MyQData>
    data_transfer(FE_Q<dim>(2), QGauss<dim>(3), rhs);
  {
    DoFHandler<dim> dof_handler(tr);
    FE_Q<dim>       dummy_fe(1);
    FEValues<dim>   fe_values(dummy_fe, rhs, update_quadrature_points);
    dof_handler.distribute_dofs(dummy_fe);
    for (cell = tr.begin_active(); cell != tr.end(); ++cell)
      if (cell->is_locally_owned())
        {
          typename DoFHandler<dim>::active_cell_iterator dof_cell(*cell,
                                                                  &dof_handler);
          fe_values.reinit(dof_cell);
          const std::vector<Point<dim>> &q_points =
            fe_values.get_quadrature_points();
          data_storage.initialize(cell, rhs.size());
          std::vector<std::shared_ptr<MyQData>> qpd =
            data_storage.get_data(cell);
          for (unsigned int q = 0; q < rhs.size(); q++)
            qpd[q]->value = my_func.value(q_points[q]);
        }
    dof_handler.clear();
  }

  check_qph(tr, data_storage, rhs, my_func);


  // mark some for refinement
  for (cell = tr.begin_active(); cell != tr.end(); ++cell)
    if (cell->center()[0] < 0.5)
      cell->set_refine_flag();

  data_transfer.prepare_for_coarsening_and_refinement(tr, data_storage);

  tr.execute_coarsening_and_refinement();

  // create qhp data
  for (cell = tr.begin_active(); cell != tr.end(); ++cell)
    if (cell->is_locally_owned())
      {
        data_storage.initialize(cell, rhs.size());
      }

  data_transfer.interpolate();

  // check that projected data still consistent with what we expect, i.e.
  // exact reproduction of the function
  check_qph(tr, data_storage, rhs, my_func);

  deallog << "Ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  test<2>();
}
