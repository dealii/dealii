// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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


// use signals to monitor solutions converging


#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

#include "../testmatrix.h"


SolverControl::State
monitor_norm(const unsigned int    iteration,
             const double          check_value,
             const Vector<double> &current_iterate)
{
  deallog << "   -- " << iteration << ' ' << check_value << std::endl;
  deallog << "   Norm=" << current_iterate.l2_norm() << std::endl;
  return SolverControl::success;
}


SolverControl::State
monitor_mean(const unsigned int    iteration,
             const double          check_value,
             const Vector<double> &current_iterate)
{
  deallog << "   Mean=" << current_iterate.mean_value() << std::endl;
  return SolverControl::success;
}



template <typename SolverType,
          typename MatrixType,
          typename VectorType,
          class PRECONDITION>
void
check_solve(SolverType &        solver,
            const MatrixType &  A,
            VectorType &        u,
            VectorType &        f,
            const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A, u, f, P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

int
main()
{
  std::ofstream logfile("output");
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);

  GrowingVectorMemory<> mem;

  SolverControl control(100, 1.e-3);

  // create CG and GMRES solvers and attach monitors to it
  SolverCG<> cg(control, mem);
  cg.connect(&monitor_norm);
  cg.connect(&monitor_mean);

  SolverGMRES<> gmres(control,
                      mem,
                      SolverGMRES<>::AdditionalData(/*max_vecs=*/8));
  gmres.connect(&monitor_norm);
  gmres.connect(&monitor_mean);

  for (unsigned int size = 4; size <= 30; size *= 3)
    {
      unsigned int dim = (size - 1) * (size - 1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix        testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double> A(structure);
      testproblem.five_point(A);

      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);

      Vector<double> f(dim);
      Vector<double> u(dim);
      Vector<double> res(dim);

      f = 1.;
      u = 1.;

      A.residual(res, u, f);
      A.SOR(res);
      res.add(1., u);
      A.SOR_step(u, f);
      res.add(-1., u);

      deallog << "SOR-diff:" << res * res << std::endl;

      try
        {
          check_solve(cg, A, u, f, prec_ssor);
          check_solve(gmres, A, u, f, prec_ssor);
        }
      catch (std::exception &e)
        {
          std::cerr << "Exception: " << e.what() << std::endl;
        }
    }
}
