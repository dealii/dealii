// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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


#include <deal.II/algorithms/operator.h>
#include <deal.II/base/logstream.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType>
  OutputOperator<VectorType>::~OutputOperator()
  {}

  template <typename VectorType>
  OutputOperator<VectorType>::OutputOperator()
    :
    os(0)
  {}

  template <typename VectorType>
  void OutputOperator<VectorType>::initialize_stream(std::ostream &stream)
  {
    os =&stream;
  }

  template <typename VectorType>
  OutputOperator<VectorType> &
  OutputOperator<VectorType>::operator<< (const AnyData &vectors)
  {
    if (os == 0)
      {
        deallog << "Step " << step << std::endl;
        for (unsigned int i=0; i<vectors.size(); ++i)
          {
            const VectorType *v = vectors.try_read_ptr<VectorType>(i);
            if (v == 0) continue;
            deallog << vectors.name(i);
            for (unsigned int j=0; j<v->size(); ++j)
              deallog << ' ' << (*v)(j);
            deallog << std::endl;
          }
        deallog << std::endl;
      }
    else
      {
        (*os) << ' ' << step;
        for (unsigned int i=0; i<vectors.size(); ++i)
          {
            const VectorType *v = vectors.try_read_ptr<VectorType>(i);
            if (v == 0) continue;
            for (unsigned int j=0; j<v->size(); ++j)
              (*os) << ' ' << (*v)(j);
          }
        (*os) << std::endl;
      }
    return *this;
  }
}

DEAL_II_NAMESPACE_CLOSE
