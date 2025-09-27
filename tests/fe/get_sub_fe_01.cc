// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test FiniteElement::get_sub_fe()

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
works(const FiniteElement<dim> &fe, const ComponentMask &m)
{
  deallog << "FE: " << fe.get_name() << " mask: " << m << std::endl;

  const FiniteElement<dim> &child = fe.get_sub_fe(m);

  deallog << "  worked: " << child.get_name() << std::endl;
}

template <int dim>
void
fails(const FiniteElement<dim> &fe, const ComponentMask &m)
{
  deallog << "FE: " << fe.get_name() << " mask: " << m << std::endl;

  try
    {
      const FiniteElement<dim> &child = fe.get_sub_fe(m);
      deallog << "  ERROR: we succeeded and got " << child.get_name()
              << " but we should have failed!" << std::endl;
    }
  catch (...)
    {
      deallog << "  failed as expected" << std::endl;
    }
}


template <int dim>
void
check()
{
  auto mask_none = [](const unsigned int n_components) -> ComponentMask {
    return ComponentMask(n_components, false);
  };
  auto mask_all = [](const unsigned int n_components) -> ComponentMask {
    return ComponentMask(n_components, true);
  };
  auto mask_single = [](const unsigned int n_components,
                        const unsigned int single) -> ComponentMask {
    ComponentMask c(n_components, false);
    c.set(single, true);
    return c;
  };
  auto mask = [](const unsigned int n_components,
                 const unsigned int first,
                 const unsigned int last) -> ComponentMask {
    ComponentMask c(n_components, false);
    for (unsigned int i = first; i <= last; ++i)
      c.set(i, true);
    return c;
  };

  FE_Q<dim>   fe_q(2);
  FE_DGP<dim> fe_dg(2);
  FE_BDM<dim> fe_nonprim(1);

  // simple FE:
  {
    works(fe_q, mask_all(1));
    fails(fe_q, mask_none(1));
  }

  // simple non-primitive:
  {
    works(fe_nonprim, ComponentMask()); // un-sized "select all"
    works(fe_nonprim, mask_all(dim));
    fails(fe_nonprim, mask_none(dim));
    fails(fe_nonprim, mask_single(dim, 0)); // can not select part of it
  }

  // remove system:
  {
    FESystem<dim> fe(fe_q, 1);
    FESystem<dim> fe2(fe, 1);
    works(fe, mask_all(1));
    works(fe2, mask_all(1));
  }

  // part of system:
  {
    FESystem<dim> fe(fe_q, 3);
    works(fe, mask_all(3));       // keep system
    works(fe, mask_single(3, 1)); // select one component
    fails(fe,
          mask(3,
               1,
               2)); // select more than one component but not the whole FESystem
  }

  // more systems:
  {
    FESystem<dim> fe(fe_nonprim, 1, fe_dg, 1);
    works(fe, mask(dim + 1, 0, dim));     // nonprimitive
    works(fe, mask_single(dim + 1, dim)); // select one component
  }
  {
    FESystem<dim> fe(fe_nonprim, 1, fe_dg, 2);
    // non-contiguous is a fail!
    auto m = mask(dim + 2, 0, dim);
    m.set(1, false);
    m.set(dim + 1, true);
    fails(fe, m);
  }
  {
    FESystem<dim> fe(fe_q, dim, fe_q, 1);
    works(fe, mask_single(dim + 1, dim));
    fails(fe, mask(dim + 1, 0, dim - 1)); // can not select the dim fe_q
  }
  {
    FESystem<dim> qs(fe_q, dim);
    FESystem<dim> fe(qs, 1, fe_q, 1);
    works(fe, mask(dim + 1, 0, dim - 1)); // but works if nested
  }

  {
    FESystem<dim> fe(fe_q, 2, fe_nonprim, 2, fe_dg, 1);
    works(fe, mask_single(2 * dim + 3, 0));
    works(fe, mask_single(2 * dim + 3, 1));
    works(fe, mask_single(2 * dim + 3, 2 * dim + 3 - 1));
    works(fe, mask(2 * dim + 3, 2, 2 + dim - 1));           // 1st nonprimitive
    works(fe, mask(2 * dim + 3, 2 + dim, 2 + 2 * dim - 1)); // 2nd nonprimitive
  }

  // nesting doll:
  {
    FESystem<dim> inner(fe_dg, 1);
    FESystem<dim> fe(fe_nonprim, 1, inner, 2);
    FESystem<dim> outer(fe, 2);
    FESystem<dim> outer_outer(outer, 1);

    works(fe, mask_single(dim + 2, dim));
    works(fe, mask_single(dim + 2, dim + 1));
    works(outer, mask_single(2 * (dim + 2), dim));
    works(outer_outer, mask_single(2 * (dim + 2), dim));
  }
}



int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  check<2>();
}
