// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


#ifndef __deal2__dof_output_operator_h
#define __deal2__dof_output_operator_h

#include <deal.II/base/config.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/named_data.h>
#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/event.h>
#include <deal.II/algorithms/operator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  /**
   * An output operator writing a separate file in each step and
   * writing the vectors as finite element functions with respect to a
   * given DoFHandler.
   */
  template <class VECTOR, int dim, int spacedim=dim>
  class DoFOutputOperator : public OutputOperator<VECTOR>
  {
  public:
    /*
     * Constructor. The <tt>filename</tt> is the common base name of
     * all files and the argument <tt>digits</tt> should be the number
     * of digits of the highest number in the sequence. File names by
     * default have the form "outputNN" with NNN the number set by the
     * last step command. Numbers with less digits are filled with
     * zeros from the left.
     */
    DoFOutputOperator (const std::string filename_base = std::string("output"),
                       const unsigned int digits = 3);

    void parse_parameters(ParameterHandler &param);
    void initialize (const DoFHandler<dim, spacedim> &dof_handler);

    virtual OutputOperator<VECTOR> &
    operator << (const AnyData &vectors);

  private:
    SmartPointer<const DoFHandler<dim, spacedim>,
                 DoFOutputOperator<VECTOR, dim, spacedim> > dof;

    const std::string filename_base;
    const unsigned int digits;

    DataOut<dim> out;
  };

  template <class VECTOR, int dim, int spacedim>
  inline void
  DoFOutputOperator<VECTOR, dim, spacedim>::initialize(const DoFHandler<dim, spacedim> &dof_handler)
  {
    dof = &dof_handler;
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
