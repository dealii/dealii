//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by Guido Kanschat
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <algorithms/operator.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR>
  Operator<VECTOR>::~Operator()
  {}

  
  
  template <class VECTOR>
  void Operator<VECTOR>::notify(const Event& e)
  {
    notifications += e;
  }
  
  
  template <class VECTOR>
  void
  Operator<VECTOR>::clear_events ()
  {
    notifications.clear();
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
  void OutputOperator<VECTOR>::initialize_stream(std::ostream& stream)
  {
    os =&stream;
  }
  
  template <class VECTOR>
  OutputOperator<VECTOR>&
  OutputOperator<VECTOR>::operator<< (const NamedData<VECTOR*>& vectors)
  {
    if(os == 0)
    {
      //TODO: make this possible
      //deallog << ' ' << step;
      //for (unsigned int i=0;i<vectors.size();++i)
      //  vectors(i)->print(deallog);
      //deallog << std::endl;
    }
    else
    {
      (*os) << ' ' << step;
      for (unsigned int i=0;i<vectors.size();++i)
	for (unsigned int j=0;j<vectors(i)->size();++j)
	  (*os) << ' ' << (*vectors(i))(j);
      (*os) << std::endl;
    }
    return *this;
  }
}

DEAL_II_NAMESPACE_CLOSE
