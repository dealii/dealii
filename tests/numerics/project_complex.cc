// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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



// Check that VectorTools::project can be instantiated for std::complex.

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_project.templates.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"

using ComplexNumber = std::complex<double>;
using Constraints   = dealii::AffineConstraints<ComplexNumber>;
using Position      = dealii::Point<3, double>;


template void
dealii::VectorTools::project<3, dealii::Vector<ComplexNumber>, 3>(
  const dealii::Mapping<3, 3> &,
  const dealii::DoFHandler<3, 3> &,
  const Constraints &,
  const dealii::Quadrature<3> &,
  const dealii::Function<3, ComplexNumber> &,
  dealii::Vector<ComplexNumber> &,
  const bool,
  const dealii::Quadrature<2> &,
  const bool);

class Field : public dealii::Function<3, ComplexNumber>
{
public:
  Field()
  {}

  ComplexNumber
  value(const Position &, const unsigned int) const
  {
    return ComplexNumber(0, 0);
  }

  void
  vector_value(const Position &, dealii::Vector<ComplexNumber> &value) const
  {
    value.reinit(3);
    for (unsigned int i = 0; i < 3; i++)
      {
        value[i] = {0, 0};
      }
  }

  dealii::Tensor<1, 3, ComplexNumber>
  curl(const Position &) const
  {
    dealii::Tensor<1, 3, ComplexNumber> ret;
    for (unsigned int i = 0; i < 3; i++)
      {
        ret[i] = {0, 0};
      }
    return ret;
  }

  dealii::Tensor<1, 3, ComplexNumber>
  val(const Position &) const
  {
    dealii::Tensor<1, 3, ComplexNumber> ret;
    for (unsigned int i = 0; i < 3; i++)
      {
        ret[i] = {0, 0};
      }
    return ret;
  }
};


int
main()
{
  initlog();

  Triangulation<3>          triangulation;
  Field                     f;
  std::vector<unsigned int> repetitions;
  repetitions.push_back(5);
  repetitions.push_back(5);
  repetitions.push_back(5);
  dealii::Point<3, double> l(0, 0, 0);
  dealii::Point<3, double> r(1, 1, 1);
  dealii::GridGenerator::subdivided_hyper_rectangle(
    triangulation, repetitions, l, r, false);

  FE_Q<3>               fe(1);
  dealii::DoFHandler<3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  Constraints local_constraints;
  local_constraints.close();

  dealii::QGauss<3>             quadrature(2);
  dealii::Vector<ComplexNumber> interpolated_exact_solution(
    dof_handler.n_dofs());
  VectorTools::project(
    dof_handler, local_constraints, quadrature, f, interpolated_exact_solution);

  deallog << "OK" << std::endl;
}
