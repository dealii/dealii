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

#include <numerics/operator.h>
#include <numerics/data_out.h>

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
      (*os) << ' ' << step << ' ';
      for (unsigned int i=0;i<vectors.size();++i)
        vectors(i)->print(*os);
      (*os) << std::endl;
    }
    return *this;
  }

  template <class VECTOR, int dim, int spacedim>
  OutputOperator<VECTOR>& 
  DoFOutputOperator<VECTOR, dim, spacedim>::operator<<(
      const NamedData<VECTOR*>& vectors)
  {
    Assert ((dof!=0), ExcNotInitialized());
    DataOut<dim> out;
    out.attach_dof_handler (*dof);
    out.add_data_vector (*vectors(vectors.find("solution")), "solution");
    out.add_data_vector (*vectors(vectors.find("update")), "update");
    out.add_data_vector (*vectors(vectors.find("residual")), "residual");
    std::ostringstream streamOut;
    streamOut << "Newton_" << std::setw(3) << std::setfill('0') << this->step;
    std::ofstream out_filename (streamOut.str().c_str());
    out.build_patches (2);
    out.write_gnuplot (out_filename);
    return *this;
  }
}

DEAL_II_NAMESPACE_CLOSE
