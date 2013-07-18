// ---------------------------------------------------------------------
// $Id$
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
#include <deal.II/numerics/data_out.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR, int dim, int spacedim>
  OutputOperator<VECTOR> &
  DoFOutputOperator<VECTOR, dim, spacedim>::operator<<(
    const NamedData<VECTOR *> &vectors)
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
