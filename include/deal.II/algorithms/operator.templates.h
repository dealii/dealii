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

#ifndef dealii_operator_templates_h
#define dealii_operator_templates_h


#include <deal.II/base/config.h>

#include <deal.II/algorithms/operator.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector_element_access.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType>
  OutputOperator<VectorType>::OutputOperator()
    : step(numbers::invalid_unsigned_int)
    , os(nullptr)
  {}

  template <typename VectorType>
  void
  OutputOperator<VectorType>::initialize_stream(std::ostream &stream)
  {
    os = &stream;
  }

  template <typename VectorType>
  OutputOperator<VectorType> &
  OutputOperator<VectorType>::operator<<(const AnyData &vectors)
  {
    if (os == nullptr)
      {
        deallog << "Step " << step << std::endl;
        for (unsigned int i = 0; i < vectors.size(); ++i)
          {
            const VectorType *v = vectors.try_read_ptr<VectorType>(i);
            if (v == nullptr)
              continue;
            deallog << vectors.name(i);
            for (unsigned int j = 0; j < v->size(); ++j)
              deallog << ' '
                      << ::dealii::internal::ElementAccess<VectorType>::get(*v,
                                                                            j);
            deallog << std::endl;
          }
        deallog << std::endl;
      }
    else
      {
        (*os) << ' ' << step;
        for (unsigned int i = 0; i < vectors.size(); ++i)
          {
            const VectorType *v = vectors.try_read_ptr<VectorType>(i);
            if (v == nullptr)
              continue;
            for (unsigned int j = 0; j < v->size(); ++j)
              (*os) << ' '
                    << ::dealii::internal::ElementAccess<VectorType>::get(*v,
                                                                          j);
          }
        (*os) << std::endl;
      }
    return *this;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
