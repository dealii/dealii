// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests DataPostprocessor: access the cell and face we are currently
// working on for a scalar finite element field and DataOutFaces. This
// is then used to output the boundary id of a face.


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



std::ofstream logfile("output");


template <int dim>
class BoundaryIds : public DataPostprocessorScalar<dim>
{
public:
  BoundaryIds()
    : DataPostprocessorScalar<dim>("boundary_id", update_quadrature_points)
  {}


  virtual void
  evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &inputs,
    std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension(computed_quantities.size(), inputs.solution_values.size());

    // Get the cell and face we are currently dealing with:
    const typename DoFHandler<dim>::active_cell_iterator cell =
      inputs.template get_cell<dim>();
    const unsigned int face = inputs.get_face_number();

    // Then fill the output fields with the boundary_id of the face
    for (auto &output : computed_quantities)
      {
        AssertDimension(output.size(), 1);
        output(0) = cell->face(face)->boundary_id();
      }
  }
};



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  FE_DGQ<dim>        fe(1);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1, true);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);

  // Create a dummy vector. We will ignore its contents.
  Vector<double> solution(dof_handler.n_dofs());

  BoundaryIds<dim>  p;
  DataOutFaces<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, p);
  data_out.build_patches();
  data_out.write_gnuplot(logfile);
}


int
main()
{
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);

  test<2>();
  test<3>();

  return 0;
}
