// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// write the data in deal.II intermediate form, read it back in, and
// make sure that the result is the same

#include "../tests.h"
#include "data_out_common.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out_rotation.h>


std::string output_file_name = "output";


// have a class that makes sure we can get at the patches and data set
// names that the base class generates
template <int dim>
class XDataOut : public DataOutRotation<dim>
{
public:
  const std::vector<typename ::DataOutBase::Patch<dim+1,dim+1> > &
  get_patches() const
  {
    return DataOutRotation<dim>::get_patches();
  }

  std::vector<std::string>
  get_dataset_names () const
  {
    return DataOutRotation<dim>::get_dataset_names();
  }
};

// have a class that makes sure we can get at the patches and data set
// names that the base class generates
template <int dim>
class XDataOutReader : public DataOutReader<dim+1>
{
public:
  const std::vector<typename ::DataOutBase::Patch<dim+1,dim+1> > &
  get_patches() const
  {
    return DataOutReader<dim+1>::get_patches();
  }

  std::vector<std::string>
  get_dataset_names () const
  {
    return DataOutReader<dim+1>::get_dataset_names();
  }
};


void
my_check_this (const DoFHandler<3> &,
               const Vector<double> &,
               const Vector<double> &)
{
  // no checks in 3d
}



template <int dim>
void
my_check_this (const DoFHandler<dim> &dof_handler,
               const Vector<double>  &v_node,
               const Vector<double>  &v_cell)
{
  XDataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (v_node, "node_data", XDataOut<dim>::type_dof_data);
  data_out.add_data_vector (v_cell, "cell_data", XDataOut<dim>::type_cell_data);
  data_out.build_patches (4);

  {
    std::ofstream tmp ("data_out_rotation_03.tmp");
    // use full precision to avoid
    // errors
    tmp.precision (20);
    data_out.write_deal_II_intermediate (tmp);
  }

  XDataOutReader<dim> reader;
  {
    std::ifstream tmp ("data_out_rotation_03.tmp");
    reader.read (tmp);
  }

  // finally make sure that we have
  // read everything back in
  // correctly
  Assert (data_out.get_dataset_names() == reader.get_dataset_names(),
          ExcInternalError());

  Assert (data_out.get_patches().size() == reader.get_patches().size(),
          ExcInternalError());

  for (unsigned int i=0; i<reader.get_patches().size(); ++i)
    Assert (data_out.get_patches()[i] == reader.get_patches()[i],
            ExcInternalError());

  // for good measure, delete tmp file
  remove ("data_out_rotation_03.tmp");

  deallog << "OK" << std::endl;
}


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node,
            const Vector<double>  &v_cell)
{
  // since we can't forward declare
  // check_this in this file (it is forward
  // declared in data_out_common.h), we
  // also can't make the driver file aware of
  // the overload for 1d. to avoid linker
  // errors, we can consequently not overload
  // check_this, and need this forwarder
  // function
  my_check_this (dof_handler, v_node, v_cell);
}
