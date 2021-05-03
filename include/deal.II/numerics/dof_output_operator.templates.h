// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_output_operator_templates_h
#define dealii_dof_output_operator_templates_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/dof_output_operator.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType, int dim, int spacedim>
  DoFOutputOperator<VectorType, dim, spacedim>::DoFOutputOperator(
    const std::string &filename_base,
    const unsigned int digits)
    : filename_base(filename_base)
    , digits(digits)
  {
    out.set_default_format(DataOutBase::gnuplot);
  }


  template <typename VectorType, int dim, int spacedim>
  void
  DoFOutputOperator<VectorType, dim, spacedim>::parse_parameters(
    ParameterHandler &param)
  {
    out.parse_parameters(param);
  }

  template <typename VectorType, int dim, int spacedim>
  OutputOperator<VectorType> &
  DoFOutputOperator<VectorType, dim, spacedim>::operator<<(const AnyData &data)
  {
    Assert((dof != nullptr), ExcNotInitialized());
    out.attach_dof_handler(*dof);
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        const VectorType *p = data.try_read_ptr<VectorType>(i);
        if (p != nullptr)
          {
            out.add_data_vector(*p, data.name(i));
          }
      }
    std::ostringstream streamOut;
    streamOut << filename_base << std::setw(digits) << std::setfill('0')
              << this->step << out.default_suffix();
    std::ofstream out_filename(streamOut.str().c_str());
    out.build_patches();
    out.write(out_filename);
    out.clear();
    return *this;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
