// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_system.h>
#include <fstream>
#include <string>


template<int dim>
void print_constant_modes(const FiniteElement<dim> &fe)
{
  deallog << "Testing " << fe.get_name() << std::endl;

  Table<2,bool> constant_modes = fe.get_constant_modes().first;
  for (unsigned int r=0; r<constant_modes.n_rows(); ++r)
    {
      for (unsigned int c=0; c<constant_modes.n_cols(); ++c)
        deallog << constant_modes(r,c) << " ";
      deallog << std::endl;
    }
  deallog << std::endl;
}


template<int dim>
void test()
{
  print_constant_modes(FE_Q<dim>(1));
  print_constant_modes(FE_Q<dim>(2));
  print_constant_modes(FE_DGQ<dim>(1));
  print_constant_modes(FE_DGP<dim>(2));
  print_constant_modes(FE_Q_Hierarchical<dim>(1));
  print_constant_modes(FE_Q_Hierarchical<dim>(2));
  print_constant_modes(FE_FaceQ<dim>(1));
  print_constant_modes(FE_FaceP<dim>(1));
  print_constant_modes(FESystem<dim>(FE_Q<dim>(1), 2, FE_Q<dim>(2), 1));
  print_constant_modes(FESystem<dim>(FE_DGP<dim>(1), 1, FE_Q_iso_Q1<dim>(2), 1));
  print_constant_modes(FE_Q_DG0<dim>(1));
  print_constant_modes(FESystem<dim>(FE_Q_DG0<dim>(2), 1, FE_Q<dim>(1), 2));
  print_constant_modes(FESystem<dim>(FE_Q<dim>(1), 2, FE_Q_DG0<dim>(1), 2));
}

template <>
void test<1>()
{
  print_constant_modes(FE_Q<1>(1));
  print_constant_modes(FESystem<1>(FE_Q<1>(1), 2, FE_Q<1>(2), 1));
  print_constant_modes(FESystem<1>(FE_DGP<1>(1), 1, FE_Q_iso_Q1<1>(2), 1));
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}



