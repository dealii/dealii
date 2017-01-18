// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// Test the ouput of MatrixFree::is_supported for various FiniteElements

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_rannacher_turek.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_p1nc.h>
#include <deal.II/fe/fe_trace.h>

#include <deal.II/fe/fe_tools.h>
#include <deal.II/matrix_free/matrix_free.h>


template <int dim, int spacedim>
void print(const FiniteElement<dim, spacedim> &fe)
{
  deallog << fe.get_name()  << " supported by MatrixFree: "
          << std::boolalpha << MatrixFree<dim>::is_supported(fe)
          << std::endl;
}


int main ()
{
  initlog();

  print<2,2> (FE_ABF<2>(0));
  print<2,2> (FE_BDM<2>(1));
  print<2,2> (FE_Bernstein<2>(1));
  print<2,2> (FE_DGBDM<2>(1));
  print<2,2> (FE_DGNedelec<2>(1));
  print<1,2> (FE_DGP<1,2>(1));
  print<2,2> (FE_DGP<2>(1));
  print<2,2> (FE_DGPMonomial<2>(1));
  print<2,2> (FE_DGPNonparametric<2>(1));
  print<2,2> (FE_DGRaviartThomas<2>(1));
  print<1,2> (FE_DGQ<1,2>(1));
  print<2,2> (FE_DGQ<2>(1));
  print<2,2> (FE_DGQArbitraryNodes<2>(QGauss<1>(2)));
  print<2,2> (FE_FaceP<2>(1));
  print<2,2> (FE_FaceQ<2>(1));
  print<2,2> (FE_P1NC());
  print<2,2> (FE_Nedelec<2>(1));
  print<2,2> (FE_Nothing<2>());
  print<1,2> (FE_Q<1,2>(1));
  print<2,2> (FE_Q<2>(1));
  print<2,2> (FE_Q_Bubbles<2>(1));
  print<2,2> (FE_Q_DG0<2>(1));
  print<2,2> (FE_Q_Hierarchical<2>(1));
  print<2,2> (FE_Q_iso_Q1<2>(1));
  print<2,2> (FE_RannacherTurek<2>(0));
  print<2,2> (FE_RaviartThomas<2>(1));
  print<2,2> (FE_RaviartThomasNodal<2>(1));
  print<2,2> (FE_TraceQ<2>(1));

  return 0;
}
