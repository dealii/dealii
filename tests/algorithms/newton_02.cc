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


#include <deal.II/algorithms/newton.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// verify that all debug vectors have the correct size

using namespace Algorithms;

template <typename VectorType, int dim>
class TestOutputOperator : public OutputOperator<VectorType>
{
public:
  void
  initialize(const DoFHandler<dim> &dof_handler)
  {
    dof = &dof_handler;
  }

  virtual OutputOperator<VectorType> &
  operator<<(const AnyData &vectors)
  {
    for (unsigned int i = 0; i < vectors.size(); ++i)
      {
        const VectorType *p = vectors.try_read_ptr<VectorType>(i);
        if (p != nullptr)
          {
            // this should be equal to dof->n_dofs() otherwise
            // DoFOutputOperator will complain
            deallog << p->size() << std::endl;
          }
      }
    return *this;
  }

private:
  SmartPointer<const DoFHandler<dim>, TestOutputOperator<VectorType, dim>> dof;
};

class ZeroResidual : public Algorithms::OperatorBase
{
public:
  virtual void
  operator()(AnyData &out, const AnyData &in)
  {
    const Vector<double> &in_vector =
      *in.entry<const Vector<double> *>("Newton iterate");
    Vector<double> &out_vector = *out.entry<Vector<double> *>(0);
    out_vector.reinit(in_vector.size());
  }
};

class IdentitySolver : public Algorithms::OperatorBase
{
public:
  virtual void
  operator()(AnyData &out, const AnyData &in)
  {
    const Vector<double> &in_vector =
      *in.entry<const Vector<double> *>("Newton residual");
    Vector<double> &out_vector = *out.entry<Vector<double> *>(0);
    out_vector                 = in_vector;
  }
};


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  TestOutputOperator<Vector<double>, dim> output_operator;
  output_operator.initialize(dofh);

  IdentitySolver solver;
  ZeroResidual   residual;

  Algorithms::Newton<Vector<double>> newton(residual, solver);
  newton.initialize(output_operator);

  AnyData in_data;
  AnyData out_data;

  Vector<double> solution(dofh.n_dofs());
  out_data.add<Vector<double> *>(&solution, "solution");

  newton.debug_vectors = true;
  newton(out_data, in_data);
}


int
main()
{
  initlog();

  test<2>();
}
