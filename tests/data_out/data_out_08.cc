// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


// This test documents two unrelated bugs in DataOut when used with a Filter (by
// deriving from DataOut):
// 1. The patch index computation in data_out.cc is wrong and causes an SIGV (or
// an Assert after adding that):
/*
466: --------------------------------------------------------
466: An error occurred in line <306> of file
</ssd/branch_port_the_testsuite/deal.II/source/numerics/data_out.cc> in function
466:     void dealii::DataOut<dim, DoFHandlerType>::build_one_patch(const
std::pair<typename dealii::DataOut_DoFData<DoFHandlerType, DoFHandlerType::
dimension, DoFHandlerType:: space_dimension>::cell_iterator, unsigned int>*,
dealii::internal::DataOut::ParallelData<DoFHandlerType:: dimension,
DoFHandlerType:: space_dimension>&, dealii::DataOutBase::Patch<DoFHandlerType::
dimension, DoFHandlerType:: space_dimension>&, dealii::DataOut<dim,
DoFHandlerType>::CurvedCellRegion,
std::vector<dealii::DataOutBase::Patch<DoFHandlerType:: dimension,
DoFHandlerType:: space_dimension> >&) [with int dim = 2, DoFHandlerType =
dealii::DoFHandler<2>, typename dealii::DataOut_DoFData<DoFHandlerType,
DoFHandlerType:: dimension, DoFHandlerType:: space_dimension>::cell_iterator =
dealii::TriaIterator<dealii::CellAccessor<2, 2> >] 466: The violated condition
was: 466:     cell_and_index->second < patches.size() 466: The name and call
sequence of the exception was: 466:     ExcInternalError()
*/
// 2. DataOut used begin_active() instead of first_cell() in two places which
// caused a wrong patch to be generated when the first active cell is not picked
// by the filter.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
class FilteredDataOut : public DataOut<dim>
{
public:
  FilteredDataOut(const unsigned int subdomain_id)
    : subdomain_id(subdomain_id)
  {}

  virtual typename DataOut<dim>::cell_iterator
  first_cell()
  {
    auto cell = this->dofs->begin_active();
    while ((cell != this->dofs->end()) &&
           (cell->subdomain_id() != subdomain_id))
      ++cell;

    return cell;
  }

  virtual typename DataOut<dim>::cell_iterator
  next_cell(const typename DataOut<dim>::cell_iterator &old_cell)
  {
    if (old_cell != this->dofs->end())
      {
        const IteratorFilters::SubdomainEqualTo predicate(subdomain_id);

        return ++(
          FilteredIterator<typename Triangulation<dim>::active_cell_iterator>(
            predicate, old_cell));
      }
    else
      return old_cell;
  }

private:
  const unsigned int subdomain_id;
};


template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);

  Vector<double> cell_data(4);
  for (unsigned int i = 0; i < 4; ++i)
    cell_data(i) = i * 1.0;

  // this should skip the first cell
  typename Triangulation<dim>::active_cell_iterator it = tria.begin_active();
  it->set_subdomain_id(1);

  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // we pick only subdomain==0 which will
  // skip the first of the four cells
  FilteredDataOut<dim> data_out(0);
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(cell_data,
                           "cell_data",
                           DataOut<dim>::type_cell_data);
  data_out.build_patches();

  data_out.write_deal_II_intermediate(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(2);

  check<2>();

  return 0;
}
