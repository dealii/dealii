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



// Check that CellDataStorage erase() and clear()


#include "../tests.h"

#ifdef DEAL_II_WITH_CXX11

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_point_data.h>

#include <fstream>

using namespace dealii;

template<int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    :
    Function<dim>(1)
  {}

  double value(const Point<dim> &p,
               const unsigned int comp=0) const
  {
    const double x = p[0];
    const double y = p[1];
    // some function we know we can project with FE_Q<dim>(2)
    return 0.5*x*x+2.1*y*y+2;
  }
};

const double default_value = 0.;

class MyQData
{
public:
  MyQData(): value(default_value) {};
  virtual ~MyQData() {};

  double value;
};

DeclException3 (ExcWrongValue,
                double, double, double,
                << arg1 <<" != "<< arg2 << " with delta = " << arg3 );


template<int dim>
void test()
{
  const MyFunction<dim> my_func;

  Triangulation<dim> tr;

  GridGenerator::subdivided_hyper_cube(tr, 2);
  tr.refine_global(1);
  typename Triangulation<dim,dim>::active_cell_iterator cell;

  // pppulate quadrature point data
  QGauss<dim> rhs(4);
  CellDataStorage<typename Triangulation<dim,dim>::cell_iterator,MyQData> data_storage;
  {
    DoFHandler<dim> dof_handler(tr);
    FE_Q<dim> dummy_fe(1);
    FEValues<dim> fe_values(dummy_fe, rhs, update_quadrature_points);
    dof_handler.distribute_dofs( dummy_fe);
    for (cell = tr.begin_active();
         cell != tr.end();
         ++cell)
      if (cell->is_locally_owned())
        {
          typename DoFHandler<dim>::active_cell_iterator dof_cell(*cell, &dof_handler);
          fe_values.reinit(dof_cell);
          const std::vector<Point<dim> > &q_points = fe_values.get_quadrature_points();
          data_storage.initialize(cell,rhs.size());
          {
            std::vector<std::shared_ptr<MyQData> > qpd = data_storage.get_data(cell);
            for (unsigned int q=0; q < rhs.size(); q++)
              qpd[q]->value = my_func.value(q_points[q]);
          }

          // do erase
          const bool erased = data_storage.erase(cell);
          Assert (erased,
                  ExcInternalError());
          // initialize with default constructor
          data_storage.initialize(cell,rhs.size());
          // check that values are now zero (see default constructor)
          {
            std::vector<std::shared_ptr<MyQData> > qpd = data_storage.get_data(cell);
            for (unsigned int q=0; q < rhs.size(); q++)
              AssertThrow(qpd[q]->value == default_value,
                          ExcWrongValue(qpd[q]->value,default_value,(qpd[q]->value-default_value)));
          }

        }
    dof_handler.clear();
  }

  // call clear
  data_storage.clear();

  deallog << "Ok" << std::endl;

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  test<2>();
}

#else
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
  deallog << "Ok" << std::endl;
}
#endif
