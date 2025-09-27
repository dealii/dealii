// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test checks whether ContinuousQuadratureDatatransfer allows
// each cell to have different number of quadrature point data.

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


using MaterialBase = TransferableQuadraturePointData;
using QuadratureStorage =
  CellDataStorage<typename Triangulation<2>::cell_iterator, MaterialBase>;

// history structs for each material
struct Mat0 : MaterialBase
{
  double x;
  virtual unsigned int
  number_of_values() const final
  {
    return 1;
  }
  virtual void
  pack_values(std::vector<double> &values) const final
  {
    AssertDimension(values.size(), number_of_values());
    values[0] = x;
  }
  virtual void
  unpack_values(const std::vector<double> &values) final
  {
    AssertDimension(values.size(), number_of_values());
    x = values[0];
  }
};

struct Mat1 : MaterialBase
{
  Point<2> pt;
  virtual unsigned int
  number_of_values() const final
  {
    return 2;
  }
  virtual void
  pack_values(std::vector<double> &values) const final
  {
    AssertDimension(values.size(), number_of_values());
    for (unsigned int d = 0; d < 2; ++d)
      values[d] = pt[d];
  }
  virtual void
  unpack_values(const std::vector<double> &values) final
  {
    AssertDimension(values.size(), number_of_values());
    for (unsigned int d = 0; d < 2; ++d)
      pt[d] = values[d];
  }
};

// This material has zero-size data
struct Mat2 : MaterialBase
{
  virtual unsigned int
  number_of_values() const final
  {
    return 0;
  }
  virtual void
  pack_values(std::vector<double> &values) const final
  {
    AssertDimension(values.size(), number_of_values());
  }
  virtual void
  unpack_values(const std::vector<double> &values) final
  {
    AssertDimension(values.size(), number_of_values());
  }
};

void
initialize_data(const Triangulation<2> &tria,
                QuadratureStorage      &storage,
                const unsigned int      n_data_points_per_cell)
{
  deallog << "Initializing quadrature cell data" << std::endl;
  for (auto cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        // It does not modify the data on already initialized cell data
        if (cell->material_id() == 0)
          storage.template initialize<Mat0>(cell, n_data_points_per_cell);
        else if (cell->material_id() == 1)
          storage.template initialize<Mat1>(cell, n_data_points_per_cell);
        else if (cell->material_id() == 2)
          storage.template initialize<Mat2>(cell, n_data_points_per_cell);
        // cells with material_id == 3 do not have an associated data structure
      }
}

void
assign_value_to_data(const Triangulation<2> &tria,
                     FEValues<2>            &fe_values,
                     QuadratureStorage      &storage,
                     const unsigned int      n_data_points_per_cell)
{
  deallog << "Assigning quadrature cell data" << std::endl;
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        if (cell->material_id() == 0)
          {
            auto data_vec = storage.template get_data<Mat0>(cell);
            // This material stores the x coord of quadrature point as data
            for (unsigned int i = 0; i < n_data_points_per_cell; ++i)
              data_vec[i]->x = fe_values.quadrature_point(i)[0];
          }
        else if (cell->material_id() == 1)
          {
            auto data_vec = storage.template get_data<Mat1>(cell);
            // This material stores the quadrature point as data
            for (unsigned int i = 0; i < n_data_points_per_cell; ++i)
              data_vec[i]->pt = fe_values.quadrature_point(i);
          }
      }
}

DeclException3(ExcWrongValue,
               double,
               double,
               double,
               << "Received " << arg1 << ", Expected " << arg2
               << ", delta = " << arg3);

void
check_data(const Triangulation<2>  &tria,
           FEValues<2>             &fe_values,
           const QuadratureStorage &storage,
           const unsigned int       n_data_points_per_cell)
{
  deallog << "Checking quadrature cell data" << std::endl;
  constexpr double eps = 1e-10;
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        if (cell->material_id() == 0)
          {
            const auto data_vec = storage.template get_data<Mat0>(cell);
            for (unsigned int i = 0; i < n_data_points_per_cell; ++i)
              {
                const double value         = data_vec[i]->x;
                const double correct_value = fe_values.quadrature_point(i)[0];
                AssertThrow(std::fabs(value - correct_value) < eps,
                            ExcWrongValue(value,
                                          correct_value,
                                          value - correct_value));
              }
          }
        else if (cell->material_id() == 1)
          {
            auto data_vec = storage.template get_data<Mat1>(cell);
            for (unsigned int i = 0; i < n_data_points_per_cell; ++i)
              {
                const double value =
                  fe_values.quadrature_point(i).distance(data_vec[i]->pt);
                const double correct_value = 0.;
                AssertThrow(std::fabs(value - correct_value) < eps,
                            ExcWrongValue(value,
                                          correct_value,
                                          value - correct_value));
              }
          }
      }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  // create a mesh with four cells
  parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);

  // assign material ids
  for (auto &cell : tria.active_cell_iterators())
    {
      const auto center = cell->center();
      if (center[0] < 0.5 && center[1] < 0.5)
        cell->recursively_set_material_id(0);
      else if (center[0] >= 0.5 && center[1] < 0.5)
        cell->recursively_set_material_id(1);
      else if (center[0] < 0.5 && center[1] >= 0.5)
        cell->recursively_set_material_id(2);
      else
        cell->recursively_set_material_id(3);
    }

  // create a history structure and populate it
  QuadratureStorage storage;
  QGauss<2>         quadrature(3);
  initialize_data(tria, storage, quadrature.size());

  // Assign some value to the cell data
  FE_Q<2>     fe(1);
  FEValues<2> fe_values(fe, quadrature, UpdateFlags::update_quadrature_points);
  assign_value_to_data(tria, fe_values, storage, quadrature.size());
  check_data(tria, fe_values, storage, quadrature.size());

  // refine all cells and transfer the data unto the new cells
  parallel::distributed::ContinuousQuadratureDataTransfer<2, MaterialBase>
    data_transfer(fe, /*mass_quadrature*/ QGauss<2>(2), quadrature);
  data_transfer.prepare_for_coarsening_and_refinement(tria, storage);
  tria.refine_global(1);
  initialize_data(tria,
                  storage,
                  quadrature.size()); // initialize newly created cells
  data_transfer.interpolate();
  deallog << "Refinement done" << std::endl;

  // Check the results after refinement
  check_data(tria, fe_values, storage, quadrature.size());

  deallog << "OK" << std::endl;

  return 0;
}
