// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


#include <deal.II/numerics/dof_output_operator.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR, int dim, int spacedim>
  DoFOutputOperator<VECTOR, dim, spacedim>::DoFOutputOperator (
    const std::string filename_base,
    const unsigned int digits)
    :
    filename_base(filename_base),
    digits(digits)
  {
    out.set_default_format(DataOutBase::gnuplot);
  }


  template <class VECTOR, int dim, int spacedim>
  void
  DoFOutputOperator<VECTOR, dim, spacedim>::parse_parameters(ParameterHandler &param)
  {
    out.parse_parameters(param);
  }

  template <class VECTOR, int dim, int spacedim>
  OutputOperator<VECTOR> &
  DoFOutputOperator<VECTOR, dim, spacedim>::operator<<(
    const AnyData &data)
  {
    Assert ((dof!=0), ExcNotInitialized());
    out.attach_dof_handler (*dof);
    for (unsigned int i=0; i<data.size(); ++i)
      {
        const VECTOR *p = data.try_read_ptr<VECTOR>(i);
        if (p!=0)
          {
            out.add_data_vector (*p, data.name(i));
          }
      }
    std::ostringstream streamOut;
    streamOut << filename_base
              << std::setw(digits) << std::setfill('0') << this->step
              << out.default_suffix();
    std::ofstream out_filename (streamOut.str().c_str());
    out.build_patches ();
    out.write (out_filename);
    out.clear ();
    return *this;
  }
}

DEAL_II_NAMESPACE_CLOSE
