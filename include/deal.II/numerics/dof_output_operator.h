// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_dof_output_operator_h
#define dealii_dof_output_operator_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/base/event.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
class ParameterHandler;
#endif

namespace Algorithms
{
  /**
   * An output operator writing a separate file in each step and writing the
   * vectors as finite element functions with respect to a given DoFHandler.
   */
  template <typename VectorType, int dim, int spacedim = dim>
  class DoFOutputOperator : public OutputOperator<VectorType>
  {
  public:
    /*
     * Constructor. The <tt>filename</tt> is the common base name of
     * all files and the argument <tt>digits</tt> should be the number
     * of digits of the highest number in the sequence. File names by
     * default have the form "outputNNN" with NNN the number set by the
     * last step command. Numbers with less digits are filled with
     * zeros from the left.
     */
    DoFOutputOperator(const std::string &filename_base = "output",
                      const unsigned int digits        = 3);

    void
    parse_parameters(ParameterHandler &param);
    void
    initialize(const DoFHandler<dim, spacedim> &dof_handler);

    virtual OutputOperator<VectorType> &
    operator<<(const AnyData &vectors) override;

  private:
    ObserverPointer<const DoFHandler<dim, spacedim>,
                    DoFOutputOperator<VectorType, dim, spacedim>>
      dof;

    const std::string  filename_base;
    const unsigned int digits;

    DataOut<dim> out;
  };

  template <typename VectorType, int dim, int spacedim>
  inline void
  DoFOutputOperator<VectorType, dim, spacedim>::initialize(
    const DoFHandler<dim, spacedim> &dof_handler)
  {
    dof = &dof_handler;
  }
} // namespace Algorithms


DEAL_II_NAMESPACE_CLOSE

#endif
