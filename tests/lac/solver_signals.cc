// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


// Connects slots to all signals in solver_cg and solver_gmres and writes all
//output to deallog.


#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>


void output_double_number(double input,const std::string &text)
{
  deallog<<text<< input<<std::endl;
}

void output_coefficients(double alpha,double beta)
{
  deallog<<"alpha: "<< alpha << " beta: "<< beta << std::endl;
}

template<class NUMBER>
void output_eigenvalues(const std::vector<NUMBER> &eigenvalues,const std::string &text)
{
  deallog<< text;
  for (unsigned int j = 0; j < eigenvalues.size(); ++j)
    {
      deallog<< ' ' << eigenvalues.at(j);
    }
  deallog << std::endl;
}

template<class NUMBER>
void output_hessenberg_matrix(const FullMatrix<NUMBER> &H,const std::string &text)
{
  deallog << text << std::endl;
  for (unsigned int i=0; i<H.m(); ++i)
    {
      for (unsigned int j=0; j<H.n(); ++j)
        deallog << H(i,j) << " ";
      deallog << std::endl;
    }
}

template<class NUMBER>
void output_arnoldi_vectors_norms(const internal::SolverGMRES::TmpVectors<Vector<NUMBER> > &tmp_vector,const std::string &text)
{
  deallog << text << std::endl;
  for (unsigned int i=0; i<tmp_vector.size(); ++i)
    deallog << tmp_vector[i].l2_norm() << std::endl;
}

template<typename SolverType, typename MatrixType, typename VectorType, class PRECONDITION>
void
check_solve(SolverType         &solver,
            const MatrixType   &A,
            VectorType         &u,
            VectorType         &f,
            const PRECONDITION &P)
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



int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  SolverControl solver_control(100, 1.e-3);

  unsigned int size=30;
  unsigned int dim = (size-1)*(size-1);

  // Make matrix
  FDMatrix testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  testproblem.five_point(A);

  Vector<double>  f(dim);
  Vector<double>  u(dim);
  f = 1.;

  try
    {
      SolverCG<> solver_cg(solver_control);
      //Attach all possible slots.
      solver_cg.connect_coefficients_slot(&output_coefficients);
      solver_cg.connect_condition_number_slot(
        std_cxx11::bind(output_double_number,std_cxx11::_1,"Condition number: "),true);
      solver_cg.connect_condition_number_slot(
        std_cxx11::bind(output_double_number,std_cxx11::_1,"Final condition number: "));
      solver_cg.connect_eigenvalues_slot(
        std_cxx11::bind(output_eigenvalues<double>,std_cxx11::_1,"Eigenvalues: "),true);
      solver_cg.connect_eigenvalues_slot(
        std_cxx11::bind(output_eigenvalues<double>,std_cxx11::_1,"Final Eigenvalues: "));
      solver_cg.solve(A,u,f,PreconditionIdentity());

      u=0;
      SolverGMRES<> solver_gmres(solver_control);
      //Attach all possible slots.
      solver_gmres.connect_condition_number_slot(
        std_cxx11::bind(output_double_number,std_cxx11::_1,"Condition number: "),true);
      solver_gmres.connect_condition_number_slot(
        std_cxx11::bind(output_double_number,std_cxx11::_1,"Final condition number: "));
      solver_gmres.connect_eigenvalues_slot(
        std_cxx11::bind(output_eigenvalues<std::complex<double> >,std_cxx11::_1,"Eigenvalues: "),true);
      solver_gmres.connect_eigenvalues_slot(
        std_cxx11::bind(output_eigenvalues<std::complex<double> >,std_cxx11::_1,"Final Eigenvalues: "));
      solver_gmres.connect_hessenberg_slot(
        std_cxx11::bind(output_hessenberg_matrix<double>,std_cxx11::_1,"Hessenberg matrix: "),false);
      solver_gmres.connect_krylov_space_slot(
        std_cxx11::bind(output_arnoldi_vectors_norms<double>,std_cxx11::_1,"Arnoldi vectors norms: "));
      solver_gmres.solve(A,u,f,PreconditionIdentity());
    }
  catch (std::exception &e)
    {
      std::cerr << "Exception: " << e.what() << std::endl;
    }

}
