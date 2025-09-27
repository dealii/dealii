// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test FiniteElement::has_support_points() for various elements,
// including FE_Nothing.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test_2d_3d(std::vector<FiniteElement<dim> *> &finite_elements)
{
  // Vector DG elements
  finite_elements.push_back(new FE_DGRaviartThomas<dim>(0));
  finite_elements.push_back(new FE_DGRaviartThomas<dim>(1));
  finite_elements.push_back(new FE_DGBDM<dim>(1));
  finite_elements.push_back(new FE_DGBDM<dim>(2));
  finite_elements.push_back(new FE_DGNedelec<dim>(0));
  finite_elements.push_back(new FE_DGNedelec<dim>(1));

  // Hdiv elements
  FE_RaviartThomas<dim> *rt0 = new FE_RaviartThomas<dim>(0);
  finite_elements.push_back(rt0);

  FE_RaviartThomas<dim> *rt1 = new FE_RaviartThomas<dim>(1);
  finite_elements.push_back(rt1);

  finite_elements.push_back(new FE_RaviartThomas<dim>(2));
  finite_elements.push_back(new FESystem<dim>(*rt1, 1, FE_DGQ<dim>(1), 1));

  finite_elements.push_back(new FE_BDM<dim>(1));
  finite_elements.push_back(new FE_BDM<dim>(2));

  // Hcurl elements
  FE_Nedelec<dim> *ned0 = new FE_Nedelec<dim>(0);
  finite_elements.push_back(ned0);
  FE_Nedelec<dim> *ned1 = new FE_Nedelec<dim>(1);
  finite_elements.push_back(ned1);
}

void test_2d_3d(std::vector<FiniteElement<1> *> /*&
                finite_elements*/)
{}

template <int dim>
void
test_finite_elements()
{
  std::vector<FiniteElement<dim> *> finite_elements;
  finite_elements.push_back(new FE_Nothing<dim>());
  finite_elements.push_back(new FE_Q<dim>(1));
  finite_elements.push_back(new FE_Q<dim>(2));
  finite_elements.push_back(new FE_Q<dim>(4));
  finite_elements.push_back(new FE_Q_Hierarchical<dim>(1));
  finite_elements.push_back(new FE_Q_Hierarchical<dim>(2));
  finite_elements.push_back(new FE_Q_Hierarchical<dim>(4));
  finite_elements.push_back(new FE_DGQ<dim>(1));
  finite_elements.push_back(new FE_DGQ<dim>(2));
  finite_elements.push_back(
    new FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapezoid<1>(), 4)));
  finite_elements.push_back(new FE_DGQ<dim>(4));
  finite_elements.push_back(new FE_DGQArbitraryNodes<dim>(QGauss<1>(3)));
  finite_elements.push_back(new FE_DGQLegendre<dim>(1));
  finite_elements.push_back(new FE_DGQLegendre<dim>(2));
  finite_elements.push_back(new FE_DGQHermite<dim>(3));
  finite_elements.push_back(new FE_DGP<dim>(1));
  finite_elements.push_back(new FE_DGP<dim>(2));
  finite_elements.push_back(new FESystem<dim>(FE_Q<dim>(2), 2));
  finite_elements.push_back(
    new FESystem<dim>(FE_Q<dim>(1), 2, FE_Q<dim>(2), 1));
  finite_elements.push_back(
    new FESystem<dim>(FE_Q<dim>(1), 2, FE_Q<dim>(2), 1, FE_Nothing<dim>(), 2));

  // Face Q elements
  finite_elements.push_back(new FE_FaceQ<dim>(0));
  finite_elements.push_back(new FE_FaceQ<dim>(1));
  finite_elements.push_back(new FE_FaceQ<dim>(3));
  // Face P elements
  finite_elements.push_back(new FE_FaceP<dim>(0));
  finite_elements.push_back(new FE_FaceP<dim>(1));
  finite_elements.push_back(new FE_FaceP<dim>(3));

  // Check vector elements in 2d and higher only
  test_2d_3d(finite_elements);

  if (dim == 2)
    {
      finite_elements.push_back(new FE_DGBDM<dim>(1));
      finite_elements.push_back(new FE_DGBDM<dim>(2));
    }
  if (dim > 1)
    {
      FE_RaviartThomasNodal<dim> *rt0 = new FE_RaviartThomasNodal<dim>(0);
      FE_RaviartThomasNodal<dim> *rt1 = new FE_RaviartThomasNodal<dim>(1);
      finite_elements.push_back(rt0);
      finite_elements.push_back(rt1);
      finite_elements.push_back(new FESystem<dim>(*rt1, 1, FE_DGQ<dim>(1), 1));
      finite_elements.push_back(
        new FESystem<dim>(*rt1, 1, FE_DGQ<dim>(1), 1, FE_Nothing<dim>(), 1));
    }

  finite_elements.push_back(new FESystem<dim>(FE_Q<dim>(3), 2));
  finite_elements.push_back(
    new FESystem<dim>(FE_Q<dim>(1), 2, FE_Q<dim>(3), 1));
  finite_elements.push_back(new FESystem<dim>(FE_Q<dim>(4), 2));

  // have systems of systems, and
  // construct hierarchies of
  // subsequently weirder elements by
  // taking each of them in turn as
  // basis of others
  finite_elements.push_back(
    new FESystem<dim>(FESystem<dim>(FE_Q<dim>(1), 2), 2));
  finite_elements.push_back(new FESystem<dim>(
    FESystem<dim>(FE_Q<dim>(1), 2), 1, FESystem<dim>(FE_DGQ<dim>(1), 2), 1));
  finite_elements.push_back(
    new FESystem<dim>(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1),
                      1,
                      FESystem<dim>(FE_Q<dim>(2), 2),
                      1,
                      FESystem<dim>(FE_DGQ<dim>(2), 2),
                      1));
  finite_elements.push_back(
    new FESystem<dim>(*finite_elements[finite_elements.size() - 3],
                      2,
                      *finite_elements[finite_elements.size() - 2],
                      1,
                      *finite_elements[finite_elements.size() - 1],
                      2));

  deallog << std::endl << "dim=" << dim << std::endl;
  for (unsigned int n = 0; n < finite_elements.size(); ++n)
    {
      FiniteElement<dim> *fe_data = finite_elements[n];
      deallog << "  " << fe_data->get_name() << std::endl;
      deallog << "    has_face_support_points="
              << (fe_data->has_face_support_points() ? "true" : "false")
              << std::endl;
    }

  // delete all FiniteElement objects
  for (unsigned int i = 0; i < finite_elements.size(); ++i)
    delete finite_elements[i];
}

int
main()
{
  initlog();

  test_finite_elements<1>();
  test_finite_elements<2>();
  test_finite_elements<3>();
}
