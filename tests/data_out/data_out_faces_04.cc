// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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



// write the data in deal.II intermediate form, read it back in, and
// make sure that the result is the same

// this is as the _03 test except that it tests our ability to read and write
// files with vector components in them

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out_faces.h>

#include "../tests.h"

#include "data_out_common.h"



// have a class that makes sure we can get at the patches and data set
// names that the base class generates
template <int dim>
class XDataOut : public DataOutFaces<dim>
{
public:
  const std::vector<typename ::DataOutBase::Patch<dim - 1, dim>> &
  get_patches() const
  {
    return DataOutFaces<dim>::get_patches();
  }

  std::vector<std::string>
  get_dataset_names() const
  {
    return DataOutFaces<dim>::get_dataset_names();
  }

  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const
  {
    return DataOutFaces<dim>::get_nonscalar_data_ranges();
  }
};

// have a class that makes sure we can get at the patches and data set
// names that the base class generates
template <int dim>
class XDataOutReader : public DataOutReader<dim - 1, dim>
{
public:
  const std::vector<typename ::DataOutBase::Patch<dim - 1, dim>> &
  get_patches() const
  {
    return DataOutReader<dim - 1, dim>::get_patches();
  }

  std::vector<std::string>
  get_dataset_names() const
  {
    return DataOutReader<dim - 1, dim>::get_dataset_names();
  }

  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const
  {
    return DataOutReader<dim - 1, dim>::get_nonscalar_data_ranges();
  }
};



void
my_check_this(const DoFHandler<1> &,
              const Vector<double> &,
              const Vector<double> &)
{
  // don't check in 1d
}



template <int dim>
void
my_check_this(const DoFHandler<dim> &dof_handler,
              const Vector<double> & v_node,
              const Vector<double> & v_cell)
{
  XDataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(v_node, "node_data", XDataOut<dim>::type_dof_data);
  data_out.add_data_vector(v_cell, "cell_data", XDataOut<dim>::type_cell_data);
  data_out.build_patches();

  {
    std::ofstream tmp("data_out_faces_04.tmp");
    data_out.write_deal_II_intermediate(tmp);
  }

  XDataOutReader<dim> reader;
  {
    std::ifstream tmp("data_out_faces_04.tmp");
    reader.read(tmp);
  }

  // finally make sure that we have
  // read everything back in
  // correctly
  AssertThrow(data_out.get_dataset_names() == reader.get_dataset_names(),
              ExcInternalError());

  AssertThrow(data_out.get_patches().size() == reader.get_patches().size(),
              ExcInternalError());

  for (unsigned int i = 0; i < reader.get_patches().size(); ++i)
    AssertThrow(data_out.get_patches()[i] == reader.get_patches()[i],
                ExcInternalError());

  deallog << data_out.get_nonscalar_data_ranges().size() << std::endl;
  Assert(data_out.get_nonscalar_data_ranges().size() ==
           reader.get_nonscalar_data_ranges().size(),
         ExcInternalError());
  for (unsigned int i = 0; i < data_out.get_nonscalar_data_ranges().size(); ++i)
    {
      deallog << std::get<0>(data_out.get_nonscalar_data_ranges()[i]) << ' '
              << std::get<1>(data_out.get_nonscalar_data_ranges()[i]) << ' '
              << std::get<2>(data_out.get_nonscalar_data_ranges()[i])
              << std::endl;
      Assert(std::get<0>(data_out.get_nonscalar_data_ranges()[i]) ==
               std::get<0>(reader.get_nonscalar_data_ranges()[i]),
             ExcInternalError());
      Assert(std::get<1>(data_out.get_nonscalar_data_ranges()[i]) ==
               std::get<1>(reader.get_nonscalar_data_ranges()[i]),
             ExcInternalError());
      Assert(std::get<2>(data_out.get_nonscalar_data_ranges()[i]) ==
               std::get<2>(reader.get_nonscalar_data_ranges()[i]),
             ExcInternalError());
    }

  // for good measure, delete tmp file
  remove("data_out_faces_04.tmp");

  deallog << "OK" << std::endl;
}



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double> & v_node,
           const Vector<double> & v_cell)
{
  // since we can't forward declare
  // check_this in this file (it is forward
  // declared in data_out_common.h), we
  // also can't make the driver file aware of
  // the overload for 1d. to avoid linker
  // errors, we can consequently not overload
  // check_this, and need this forwarder
  // function
  my_check_this(dof_handler, v_node, v_cell);
}
