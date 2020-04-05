// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// Check that interpolate_based_on_material_id() really does ignore
// material_ids not included in the provided std::map

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



namespace PhaseField
{
  typedef Vector<double> vectorType;


  ///////////////////////////////////////////////////////////////

  template <int dim>
  class Solid
  {
  public:
    Solid();

    ~Solid();

    void
    run();

  private:
    void
    make_grid();

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;
  };


  // constructor
  template <int dim>
  Solid<dim>::Solid()
    : fe(1)
    , dof_handler(triangulation)
  {}

  // destructor
  template <int dim>
  Solid<dim>::~Solid()
  {
    dof_handler.clear();
  }


  ///////////////////////////////////////////////////////
  template <int dim>
  void
  Solid<dim>::make_grid()
  {
    std::vector<unsigned int> repetitions(dim, 4);
    if (dim == 3)
      repetitions[dim - 2] = 4;

    GridGenerator::subdivided_hyper_rectangle(
      triangulation,
      repetitions,
      Point<dim>(),
      (dim == 1 ?
         Point<dim>(1.0) :
         (dim == 2 ? Point<dim>(1.0, 1.0) : Point<dim>(1.0, 1.0, 1.0))),
      true);


    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)

      if (cell->center()[0] < 0.5)
        cell->set_material_id(1);
      else
        cell->set_material_id(0);
  }

  //////////////////////////////////////////////////
  template <int dim>
  void
  Solid<dim>::run()
  {
    make_grid();
    dof_handler.distribute_dofs(fe);

    vectorType dst(dof_handler.n_dofs());

    deallog << " Initial values: " << std::endl;
    dst.print(deallog.get_file_stream());

    std::map<types::material_id, const Function<dim> *> function_map;

    const Functions::ConstantFunction<dim> constant_function_1(1.0);
    function_map[1] = &constant_function_1;

    VectorTools::interpolate_based_on_material_id(MappingQGeneric<dim>(1),
                                                  dof_handler,
                                                  function_map,
                                                  dst);


    deallog << " Final values: " << std::endl;
    dst.print(deallog.get_file_stream());
  }
} // namespace PhaseField


int
main()
{
  initlog();

  deallog.push("1d");
  PhaseField::Solid<1>().run();
  deallog.pop();
  deallog.push("2d");
  PhaseField::Solid<2>().run();
  deallog.pop();
  deallog.push("3d");
  PhaseField::Solid<3>().run();
  deallog.pop();
}
