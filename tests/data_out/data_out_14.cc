// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test checks whether DataOut can be copied/moved and used safely and
// no dangling references/pointers to the copied-from object remain.
// DataOut is created in functions to ensure the original object is destroyed.
// DataOut is wrapped to ensure it is copied on return, otherwise copies might
// be elided through Named Return Value Optimization (NRVO) or
// prvalue semantics.
// see https://en.cppreference.com/w/cpp/language/copy_elision.html.

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "data_out_common.h"

namespace
{
  template <int dim>
  struct Wrapper
  {
    const DataOut<dim> &
    unwrap() const
    {
      return data_out;
    }
    DataOut<dim> data_out;
  };

  template <int dim>
  Wrapper<dim>
  make_data_out()
  {
    DataOut<dim> data_out;
    // DataOut copy, not elided since it is not returned directly:
    return Wrapper<dim>{data_out};
  }

  template <int dim>
  Wrapper<dim>
  make_data_out(const DoFHandler<dim> &dof_handler)
  {
    DataOut<dim> data_out =
      make_data_out<dim>().unwrap(); // DataOut copy construction.
    data_out.attach_dof_handler(dof_handler);
    // DataOut copy, not elided since it is not returned directly:
    return Wrapper<dim>{data_out};
  }

  template <int dim>
  Wrapper<dim>
  make_data_out(const DoFHandler<dim> &dof_handler,
                const Vector<double>  &v_node,
                const Vector<double>  &v_cell)
  {
    DataOut<dim> data_out =
      make_data_out<dim>(dof_handler).unwrap(); // DataOut copy construction.
    data_out.add_data_vector(v_node, "node_data", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(v_cell, "cell_data", DataOut<dim>::type_cell_data);
    // DataOut move, make sure moves are also safe:
    return Wrapper<dim>{std::move(data_out)};
  }
} // namespace

template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double>  &v_node,
           const Vector<double>  &v_cell)
{
  DataOut<dim> data_out = make_data_out<dim>(dof_handler, v_node, v_cell)
                            .unwrap(); // DataOut copy construction.
  data_out.build_patches();
  data_out.write_dx(deallog.get_file_stream());
}
