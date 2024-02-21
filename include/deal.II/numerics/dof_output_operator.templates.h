// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_output_operator_templates_h
#define dealii_dof_output_operator_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/parameter_handler.h>

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
    std::ofstream out_filename(streamOut.str());
    out.build_patches();
    out.write(out_filename);
    out.clear();
    return *this;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
