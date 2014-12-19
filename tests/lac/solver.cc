// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/precondition.h>

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_solve( SOLVER &solver, const MATRIX &A,
             VECTOR &u, VECTOR &f, const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A,u,f,P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_Tsolve(SOLVER &solver, const MATRIX &A,
             VECTOR &u, VECTOR &f, const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.Tsolve(A,u,f,P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

int main()
{
  std::ofstream logfile("output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  GrowingVectorMemory<> mem;
  SolverControl control(100, 1.e-3);
  SolverControl verbose_control(100, 1.e-3, true);
  SolverCG<> cg(control, mem);
  SolverGMRES<> gmres(control, mem, 8);
  SolverGMRES<>::AdditionalData data(8, true);
  SolverGMRES<> gmresright(control, mem, data);
  SolverMinRes<> minres(control, mem);
  SolverBicgstab<> bicgstab(control, mem);
  SolverRichardson<> rich(control, mem);
  SolverQMRS<> qmrs(control, mem);

  for (unsigned int size=4; size <= 30; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.five_point(A);

      PreconditionIdentity prec_no;
      PreconditionRichardson prec_richardson;
      prec_richardson.initialize(0.6);
      PreconditionSOR<> prec_sor;
      prec_sor.initialize(A, 1.2);
      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);

      std::vector<types::global_dof_index> permutation(dim);
      std::vector<types::global_dof_index> inverse_permutation(dim);

      // Create a permutation: Blocks
      // backwards and every second
      // block backwards
      unsigned int k = 0;
      for (unsigned int i=0; i<size-1; ++i)
        for (unsigned int j=0; j<size-1; ++j)
          {
            if (i % 2)
              permutation[k++] = (size-i-2) * (size-1) + j;
            else
              permutation[k++] = (size-i-2) * (size-1) + size-j-2;
          }


      for (unsigned int i=0; i<dim; ++i)
        inverse_permutation[permutation[i]] = i;

      PreconditionPSOR<> prec_psor;
      prec_psor.initialize(A, permutation, inverse_permutation, 1.2);

      Vector<double>  f(dim);
      Vector<double>  u(dim);
      Vector<double> res(dim);

      f = 1.;
      u = 1.;

      A.residual(res,u,f);
      A.SOR(res);
      res.add(1.,u);
      A.SOR_step(u,f);
      res.add(-1.,u);

      deallog << "SOR-diff:" << res *res << std::endl;

      try
        {
          deallog.push("no-fail");

          control.set_max_steps(10);
          check_solve(cg,A,u,f,prec_no);
          check_solve(bicgstab,A,u,f,prec_no);
          check_solve(gmres,A,u,f,prec_no);
          check_solve(gmresright,A,u,f,prec_no);
//    check_solve(minres,A,u,f,prec_no);
          check_solve(qmrs,A,u,f,prec_no);
          control.set_max_steps(100);

          deallog.pop();

          deallog.push("no");

          rich.set_omega(1./A.diag_element(0));
          check_solve(rich,A,u,f,prec_no);
          check_solve(cg,A,u,f,prec_no);
          check_solve(bicgstab,A,u,f,prec_no);
          check_solve(gmres,A,u,f,prec_no);
          check_solve(gmresright,A,u,f,prec_no);
          check_solve(qmrs,A,u,f,prec_no);
          rich.set_omega(1.);

          deallog.pop();

          deallog.push("rich");

          rich.set_omega(1./A.diag_element(0));
          check_solve(rich,A,u,f,prec_richardson);
          check_solve(cg,A,u,f,prec_richardson);
          check_solve(bicgstab,A,u,f,prec_richardson);
          check_solve(gmres,A,u,f,prec_richardson);
          check_solve(gmresright,A,u,f,prec_richardson);
          check_solve(qmrs,A,u,f,prec_richardson);
          rich.set_omega(1.);

          deallog.pop();

          deallog.push("ssor");

          check_Tsolve(rich,A,u,f,prec_ssor);
          check_solve(rich,A,u,f,prec_ssor);
          check_solve(cg,A,u,f,prec_ssor);
          check_solve(bicgstab,A,u,f,prec_ssor);
          check_solve(gmres,A,u,f,prec_ssor);
          check_solve(gmresright,A,u,f,prec_ssor);
          check_solve(qmrs,A,u,f,prec_ssor);

          deallog.pop();

          deallog.push("sor");

          check_Tsolve(rich,A,u,f,prec_sor);
          check_solve(rich,A,u,f,prec_sor);
          check_solve(cg,A,u,f,prec_sor);
          check_solve(bicgstab,A,u,f,prec_sor);
          check_solve(gmres,A,u,f,prec_sor);
          check_solve(gmresright,A,u,f,prec_sor);

          deallog.pop();

          deallog.push("psor");

          check_Tsolve(rich,A,u,f,prec_psor);
          check_solve(rich,A,u,f,prec_psor);
          check_solve(cg,A,u,f,prec_psor);
          check_solve(bicgstab,A,u,f,prec_psor);
          check_solve(gmres,A,u,f,prec_psor);
          check_solve(gmresright,A,u,f,prec_psor);

          deallog.pop();
        }
      catch (std::exception &e)
        {
          std::cerr << "Exception: " << e.what() << std::endl;
        }
    };

  // Solve advection problem
  for (unsigned int size=4; size <= 3; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.upwind(A, true);

      PreconditionSOR<> prec_sor;
      prec_sor.initialize(A, 1.);

      std::vector<types::global_dof_index> permutation(dim);
      std::vector<types::global_dof_index> inverse_permutation(dim);

      // Create a permutation: Blocks
      // backwards and every second
      // block backwards
      unsigned int k = 0;
      for (unsigned int i=0; i<size-1; ++i)
        for (unsigned int j=0; j<size-1; ++j)
          {
            permutation[k++] = i * (size-1) + size-j-2;
          }

      for (unsigned int i=0; i<permutation.size(); ++i)
        std::cerr << ' ' << permutation[i];
      std::cerr << std::endl;

      for (unsigned int i=0; i<permutation.size(); ++i)
        inverse_permutation[permutation[i]] = i;

      for (unsigned int i=0; i<permutation.size(); ++i)
        std::cerr << ' ' << inverse_permutation[i];
      std::cerr << std::endl;

      PreconditionPSOR<> prec_psor;
      prec_psor.initialize(A, permutation, inverse_permutation, 1.);

      Vector<double>  f(dim);
      Vector<double>  u(dim);
      f = 1.;
      u = 1.;

      std::cerr << "******************************" << std::endl;

      check_solve(rich,A,u,f,prec_sor);
      check_solve(rich,A,u,f,prec_psor);
    }
}

