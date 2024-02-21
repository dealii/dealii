// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check equivalence of get_unit_support_points() and
// unit_support_point(), as well as the case where an element does not
// have support points

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_trace.h>

#include "../tests.h"


template <int dim>
void
check_cell1(const FiniteElement<dim> &fe)
{
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    AssertThrow(fe.get_unit_support_points()[i] == fe.unit_support_point(i),
                ExcInternalError());
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void
check_face1(const FiniteElement<dim> &fe)
{
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    AssertThrow(fe.get_unit_face_support_points()[i] ==
                  fe.unit_face_support_point(i),
                ExcInternalError());
  deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void
check_face1(const FiniteElement<1> &)
{}



template <int dim>
void
check1(const FiniteElement<dim> &fe)
{
  check_cell1(fe);
  check_face1(fe);
}


template <int dim>
void
check_cell2(const FiniteElement<dim> &fe, const unsigned int comp)
{
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    if (fe.system_to_component_index(i).first == comp)
      deallog << i << ' ' << fe.unit_support_point(i) << std::endl;
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void
check_face2(const FiniteElement<dim> &fe, const unsigned int comp)
{
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    if (fe.face_system_to_component_index(i).first == comp)
      deallog << i << ' ' << fe.unit_face_support_point(i) << std::endl;
  deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void
check_face2(const FiniteElement<1> &, const unsigned int)
{}



template <int dim>
void
check2(const FiniteElement<dim> &fe, const unsigned int comp)
{
  check_cell2(fe, comp);
  check_face2(fe, comp);
}



template <int dim>
void
face_check()
{
  check1(FE_FaceQ<dim>(2));
  check2(FE_FaceQ<dim>(2), 0);
  check1(FE_TraceQ<dim>(2));
  check2(FE_TraceQ<dim>(2), 0);
}



template <>
void
face_check<1>()
{}



template <int dim>
void
check()
{
  check1(FE_Q<dim>(2));
  check1(FE_DGQ<dim>(2));
  check1(FESystem<dim>(FE_Q<dim>(2), 2, FE_DGQ<dim>(2), 1));

  check1(FESystem<dim>(FE_Q<dim>(2), 2, FE_DGQ<dim>(2), 1));

  check2(FESystem<dim>(FE_Q<dim>(2), 1, FE_DGQ<dim>(2), 1), 1);
  check2(FESystem<dim>(FE_Q<dim>(2), 1, FE_DGP<dim>(2), 1), 0);
  face_check<dim>();
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
  return 0;
}
