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



// Check CellDataStorage initialize(cell,number), initialize(start,end,number),
// get_data(), and try_get_data() functions.


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

class MyQData
{
public:
  MyQData(){};
  virtual ~MyQData(){};

  double value;
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
check_qph(Triangulation<dim> &tr,
          CellDataStorage<typename Triangulation<dim, dim>::cell_iterator, DATA>
            &                    manager,
          const Quadrature<dim> &rhs_quadrature,
          const MyFunction<dim> &func)
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
        const auto &                             manager_const = manager;
        const std::vector<std::shared_ptr<DATA>> qpd = manager.get_data(cell);
        const std::vector<std::shared_ptr<const DATA>> qpd_const =
          manager_const.get_data(cell);
        const std_cxx17::optional<std::vector<std::shared_ptr<DATA>>>
          qpd_optional = manager.try_get_data(cell);
        const std_cxx17::optional<std::vector<std::shared_ptr<const DATA>>>
          qpd_const_optional = manager_const.try_get_data(cell);
        AssertThrow(qpd_optional, ExcInternalError());
        AssertThrow(qpd_const_optional, ExcInternalError());
        for (unsigned int q = 0; q < q_points.size(); q++)
          {
            const double correct_value = func.value(q_points[q]);
            const double value         = qpd[q]->value;
            AssertThrow(std::fabs(correct_value - value) < eps,
                        ExcWrongValue(correct_value,
                                      value,
                                      correct_value - value));
            const double value_const = qpd_const[q]->value;
            AssertThrow(std::fabs(correct_value - value_const) < eps,
                        ExcWrongValue(correct_value,
                                      value_const,
                                      correct_value - value_const));
            const double value_optional = (*qpd_optional)[q]->value;
            AssertThrow(std::fabs(correct_value - value_optional) < eps,
                        ExcWrongValue(correct_value,
                                      value_optional,
                                      correct_value - value_optional));
            const double value_const_optional = (*qpd_const_optional)[q]->value;
            AssertThrow(std::fabs(correct_value - value_const_optional) < eps,
                        ExcWrongValue(correct_value,
                                      value_optional,
                                      correct_value - value_const_optional));
          }
      }
  dof_handler.clear();
}

template <int dim>
void
test()
{
  const MyFunction<dim> my_func;
  Triangulation<dim>    tr;

  GridGenerator::subdivided_hyper_cube(tr, 2);
  tr.refine_global(1);
  typename Triangulation<dim, dim>::active_cell_iterator cell;

  // pppulate quadrature point data
  QGauss<dim> rhs(4);
  CellDataStorage<typename Triangulation<dim, dim>::cell_iterator, MyQData>
    data_storage;
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
          {
            // before initialization of the next cell, try_get_data must
            // return null. This test has to be done after initialize()
            // has been called at least one time
            auto next_cell = cell;
            next_cell++;
            if (next_cell != tr.end())
              {
                const std_cxx17::optional<std::vector<std::shared_ptr<MyQData>>>
                  nonexisting_data = data_storage.try_get_data(next_cell);
                AssertThrow(!nonexisting_data, ExcInternalError());
              }
          }
        }
    dof_handler.clear();
  }

  check_qph(tr, data_storage, rhs, my_func);

  data_storage.initialize(tr.begin_active(), tr.end(), rhs.size());

  check_qph(tr, data_storage, rhs, my_func);

  deallog << "Ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  test<2>();
}
