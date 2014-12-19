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


#include <deal.II/algorithms/operator.h>
#include <deal.II/base/logstream.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR>
  Operator<VECTOR>::Operator()
    : silent_compatibility(false), compatibility_flag(false)
  {}


  template <class VECTOR>
  void
  Operator<VECTOR>::operator() (AnyData &out, const AnyData &in)
  {
    // Had this function been overloaded in a derived clas, it would
    // not have been called. Therefore, we have to start the
    // compatibility engine. But before, we have to avoid an endless loop.
    Assert(!compatibility_flag, ExcMessage("Compatibility resolution of Operator generates and endless loop\n"
                                           "Please provide an operator() in a derived class"));
    compatibility_flag = true;

    NamedData<VECTOR *> new_out;
    for (unsigned int i=0; i<out.size(); ++i)
      {
        if (out.is_type<VECTOR *>(i))
          new_out.add(out.entry<VECTOR *>(i), out.name(i));
        else if (!silent_compatibility)
          deallog << "Cannot convert AnyData argument " << out.name(i) << " to NamedData"
                  << std::endl;
      }

    NamedData<VECTOR *> new_in;
    for (unsigned int i=0; i<in.size(); ++i)
      {
        //  deallog << "Convert " << in.name(i) << std::endl;
        if (in.is_type<VECTOR *>(i))
          {
            // This const cast is due to the wrong constness handling
            // in NamedData. As soon as NamedData is gone, this code
            // will not be necessary anymore. And deprecating begins
            // now.
            VECTOR *p = const_cast<VECTOR *> (in.entry<VECTOR *>(i));
            new_in.add(p, in.name(i));
          }
        else if (in.is_type<const VECTOR *>(i))
          {
            // This const cast is due to the wrong constness handling
            // in NamedData. As soon as NamedData is gone, this code
            // will not be necessary anymore. And deprecating begins
            // now.
            VECTOR *p = const_cast<VECTOR *> (in.entry<const VECTOR *>(i));
            new_in.add(p, in.name(i));
          }
        else if (!silent_compatibility)
          deallog << "Cannot convert AnyData argument " << in.name(i)
                  << " to NamedData" << std::endl;
      }
    this->operator() (new_out, new_in);
    compatibility_flag = false;
  }


  template <class VECTOR>
  void
  Operator<VECTOR>::operator() (NamedData<VECTOR *> &out, const NamedData<VECTOR *> &in)
  {
    // Had this function been overloaded in a derived clas, it would
    // not have been called. Therefore, we have to start the
    // compatibility engine. But before, we have to avoid an endless loop.
    Assert(!compatibility_flag, ExcMessage("Compatibility resolution of Operator generates and endless loop\n"
                                           "Please provide an operator() in a derived class"));
    compatibility_flag = true;

    AnyData new_out;
    for (unsigned int i=0; i<out.size(); ++i)
      new_out.add(out(i), out.name(i));

    AnyData new_in;
    for (unsigned int i=0; i<in.size(); ++i)
      new_in.add(in(i), in.name(i));

    this->operator() (new_out, new_in);
    compatibility_flag = false;
  }


  template <class VECTOR>
  OutputOperator<VECTOR>::~OutputOperator()
  {}

  template <class VECTOR>
  OutputOperator<VECTOR>::OutputOperator()
    :
    os(0)
  {}

  template <class VECTOR>
  void OutputOperator<VECTOR>::initialize_stream(std::ostream &stream)
  {
    os =&stream;
  }

  template <class VECTOR>
  OutputOperator<VECTOR> &
  OutputOperator<VECTOR>::operator<< (const AnyData &vectors)
  {
    if (os == 0)
      {
        deallog << "Step " << step << std::endl;
        for (unsigned int i=0; i<vectors.size(); ++i)
          {
            const VECTOR *v = vectors.try_read_ptr<VECTOR>(i);
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
            const VECTOR *v = vectors.try_read_ptr<VECTOR>(i);
            if (v == 0) continue;
            for (unsigned int j=0; j<v->size(); ++j)
              (*os) << ' ' << (*v)(j);
          }
        (*os) << std::endl;
      }
    return *this;
  }


  template <class VECTOR>
  OutputOperator<VECTOR> &
  OutputOperator<VECTOR>::operator<< (const NamedData<VECTOR *> &vectors)
  {
    const AnyData newdata = vectors;
    (*this) << newdata;
    return *this;
  }
}

DEAL_II_NAMESPACE_CLOSE
