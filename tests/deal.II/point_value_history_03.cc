// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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




// ************************************************
// A test program for the PointValueHistory class
// Currently this only tests a finite element system
// with 3 components, on a hyper cube and with Vectors.
// The testing is done in two dimensions.
// ************************************************



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>

#include <deal.II/numerics/point_value_history.h>

#include <fstream>

using namespace dealii;


template <int dim>
class Postprocess : public DataPostprocessor<dim>
{
public:

  void compute_derived_quantities_scalar (
    const std::vector< double > &,
    const std::vector< Tensor< 1, dim > > &,
    const std::vector< Tensor< 2, dim > > &,
    const std::vector< Point< dim > > &,
    const std::vector< Point< dim > > &,
    std::vector< Vector< double > > &
  ) const;

  std::vector<std::string> get_names () const;
  UpdateFlags              get_needed_update_flags () const;
  unsigned int             n_output_variables () const;
  // The following function is not required
  // by the point_value_history class.
  //std::vector<DataComponentInterpretation::DataComponentInterpretation>
  //                  get_data_component_interpretation () const;
};

template <int dim>
std::vector<std::string>
Postprocess<dim>::get_names() const
{
  std::vector<std::string> names;
  names.push_back ("X_gradient");
  names.push_back ("X_hessian");
  return names;
}

template <int dim>
UpdateFlags
Postprocess<dim>::get_needed_update_flags () const
{
  return update_values | update_gradients | update_hessians;
}

template <int dim>
unsigned int
Postprocess<dim>::n_output_variables () const
{
  return 2;
}


template <int dim>
void
Postprocess<dim>::compute_derived_quantities_scalar (
  const std::vector< double >                  &uh,
  const std::vector< Tensor< 1, dim > >   &duh,
  const std::vector< Tensor< 2, dim > >   &dduh,
  const std::vector< Point< dim > >                     & /* normals */,
  const std::vector< Point< dim > >                     & /* locations */,
  std::vector< Vector< double > >                        &computed_quantities
) const
{
  Assert(computed_quantities.size() == uh.size(),
         ExcDimensionMismatch (computed_quantities.size(), uh.size()));

  for (unsigned int i=0; i<computed_quantities.size(); i++)
    {
      Assert(computed_quantities[i].size() == 2,
             ExcDimensionMismatch (computed_quantities[i].size(), 2));

      computed_quantities[i](0) = duh[i][0]; // norm of x gradient
      computed_quantities[i](1) = dduh[i][0].norm(); // norm of x hessian
    }
}



template <int dim>
class TestPointValueHistory
{
public:
  TestPointValueHistory();
  void run();

private:
  void output_results (unsigned int step, Vector <double> solution) const;

  Triangulation <dim> triangulation;
  FE_Q<dim> finite_element;
  DoFHandler <dim > dof_handler;
  PointValueHistory <dim> test_copy;
  std::vector <Point <dim> > postprocessor_locations;
};




template <int dim>
TestPointValueHistory<dim>::TestPointValueHistory() :
  finite_element(2),
  dof_handler(triangulation)
{ }




template <int dim>
void TestPointValueHistory<dim>::run()
{
  // Make a triangulation
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(2); // refine 2 times to make 5 nodes per side

  // make a DOF handler, a model solution filled with ones and a flow vector
//    DoFHandler <dim > dof_handler(triangulation);
  dof_handler.distribute_dofs(finite_element);
  DoFRenumbering::Cuthill_McKee(dof_handler);

  // Vector
  Vector <double> solution, post_processed, poles;
  solution.reinit(dof_handler.n_dofs());
  post_processed.reinit(dof_handler.n_dofs());
  poles.reinit(dof_handler.n_dofs());

//            // BlockVector
  // BlockVector not used with scalar fields generally

  // set up a simple linear discrete time system so that time plots vary
  // over the mesh but can be easily checked. The basic idea is to have each
  // component of the fe_system to depend on a specific dimension (i.e component 0
  // depends on dim 0, etc. % dim handles the case where there are more components
  // than dimensions. The code breaks down at the edges of the mesh and this is
  // not corrected for. The code used in this test is simplified from point_value_history_01.
  {
    Quadrature<dim> quadrature_formula(finite_element.get_unit_support_points ());
    FEValues<dim> fe_values(finite_element, quadrature_formula, update_values | update_quadrature_points); // just need local_dof_indices and quadrature_points

    std::vector<types::global_dof_index> local_dof_indices(finite_element.dofs_per_cell);
    std::vector<Point<dim> > dof_locations(finite_element.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell, endc;
    cell = dof_handler.begin_active();
    endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell); // need to get local_dof_indices
        cell->get_dof_indices(local_dof_indices);
        dof_locations = fe_values.get_quadrature_points();

        for (unsigned int dof = 0; dof != finite_element.dofs_per_cell; dof++)
          {
            unsigned int dof_component = finite_element.system_to_component_index(dof).first;
            // The line above will always evaluate to 1,
            // simplifying the remaining lines

            poles(local_dof_indices[dof]) = -dof_locations[dof](dof_component % dim);

            if (dof_component == dim) // components start numbering at 0
              poles(local_dof_indices[dof]) = -0.1; // dim+1th component is not handled well by the code above

            solution(local_dof_indices [dof]) = 1; // start all solutions at 1
          }
      } // loop over all cells

    poles.add(1.0); // slow down the pole settling time
    post_processed.add(3.0); // set to starting value.
//        output_results(10, poles);
  }

  // Setup monitor node to print variation over time
  unsigned int n_inputs = 1;
  PointValueHistory<dim> node_monitor (dof_handler, n_inputs);
  PointValueHistory<dim> no_dof_handler (n_inputs);

  // check that the assignment operator is valid
  test_copy = node_monitor;
  test_copy.add_point(Point<2>(1, 0.2));
  test_copy.add_field_name("Solution");
  std::vector < std::vector <Point <dim> > > selected_locations;
  test_copy.get_support_locations(selected_locations);
  test_copy.mark_support_locations();
  test_copy.close();
  test_copy.start_new_dataset(0.1);
  test_copy.evaluate_field("Solution", solution);
  std::vector <double> input_value(n_inputs, 1);
  test_copy.push_back_independent(input_value);
  test_copy.write_gnuplot("Test_Copy");
  test_copy.status(deallog.get_file_stream());
  test_copy.clear ();
  // end of assignment operator check

  {
    node_monitor.add_field_name("Solution");
    std::vector <std::string> solution_names;
    solution_names.push_back("Solution");
    node_monitor.add_component_names ("Solution", solution_names);
    node_monitor.add_field_name("Post Processed Vector"); // not sensitive to spaces
    node_monitor.add_field_name("Pressure", 1);
    std::vector <bool> component_mask (1, true);
    node_monitor.add_field_name("Req_sol", component_mask);
    component_mask = std::vector <bool> (2, false);
    component_mask[0] = true;
    node_monitor.add_field_name("X_gradient", component_mask);
    component_mask = std::vector <bool> (2, false);
    component_mask[1] = true;
    node_monitor.add_field_name("X_hessian", component_mask);
    std::vector <std::string> indep_names;
    indep_names.push_back ("Input");
    node_monitor.add_independent_names(indep_names);

    // two alternatives here, adding a point at a time or a vector of points
    // 2d points
    std::vector <Point < 2 > > point_vector(5, Point < 2 > ());
    point_vector[0] = Point < 2 > (0, 0); // some of these points will hit a node, others won't
    point_vector[1] = Point < 2 > (0.25, 0);
    point_vector[2] = Point < 2 > (0.25, 0.45);
    point_vector[3] = Point < 2 > (0.45, 0.45);
    point_vector[4] = Point < 2 > (0.8, 0.8);

    node_monitor.add_points(point_vector);
    node_monitor.add_point(Point<2>(1, 0.2)); // add a single point

    // MonitorNode requires that the instance is 'closed' before any data is added
    // this ensures that points are not added once time starts.
    node_monitor.close();
    no_dof_handler.close(); // closing still required!

    std::vector < std::vector <Point <dim> > > selected_locations;
    node_monitor.get_support_locations(selected_locations);
    Vector<double> node_locations = node_monitor.mark_support_locations();
    QGauss<dim> postprocess_quadrature (2);
    node_monitor.get_postprocessor_locations(postprocess_quadrature, postprocessor_locations);
  }

  double delta_t = 0.000001;
  double t_max = 0.00001;
  unsigned int step = 0;

  for (double time = 0; time < t_max; time = time + delta_t)
    {
      node_monitor.start_new_dataset(time);
      no_dof_handler.start_new_dataset(time);
      // time and input are special, they don't vary over the mesh.

      std::vector <double> input_value(n_inputs, 1); // manufacture an input value
      node_monitor.push_back_independent(input_value);
      no_dof_handler.push_back_independent(input_value);
      node_monitor.evaluate_field("Solution", solution);
      node_monitor.evaluate_field("Post Processed Vector", post_processed);
      node_monitor.evaluate_field("Pressure", solution);
      node_monitor.evaluate_field_at_requested_location("Req_sol", solution);

      Postprocess<dim> postprocessor;
      QGauss<dim> postprocess_quadrature (2);
      std::vector<std::string> names;
      names.push_back ("X_gradient");
      names.push_back ("X_hessian");
      node_monitor.evaluate_field(names, solution, postprocessor, postprocess_quadrature);

//        output_results (step, solution);
      step++;

      solution.scale(poles); // decaying exponentials of varying time constants
      post_processed = solution;
      post_processed.add(2.0); // simple post processing, giving it a dc offset
    }
  node_monitor.write_gnuplot("node", postprocessor_locations);
  no_dof_handler.write_gnuplot("no_dof"); // no point in adding postprocessor_locations

  node_monitor.status (deallog.get_file_stream());
  no_dof_handler.status (deallog.get_file_stream());

  deallog << "Starting data files" << std::endl;

  // copy all the data into deallog and
  // delete those files
  const std::string filenames[]
    = { "node_00.gpl",
        "node_01.gpl",
        "node_02.gpl",
        "node_03.gpl",
        "node_04.gpl",
        "node_05.gpl",
        "node_indep.gpl",
        "Test_Copy_00.gpl",
        "Test_Copy_indep.gpl",
        "no_dof_indep.gpl"
      };

  for (unsigned int i=0; i<sizeof(filenames)/sizeof(filenames[0]); ++i)
    {
      deallog << "Copying output file " << filenames[i]
              << std::endl;

      std::ifstream in(filenames[i].c_str());
      AssertThrow (in, ExcIO());

      std::string s;
      while (in)
        {
          std::getline (in, s);
          deallog << s << std::endl;
        }

      std::remove (filenames[i].c_str());

      deallog << std::endl;
    }
}

template <int dim>
void TestPointValueHistory<dim>::output_results (unsigned int step, Vector <double> solution)  const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches (2);

  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (step, 2)
           << ".gpl";

  std::ofstream output (filename.str().c_str());
  data_out.write_gnuplot (output);
}



int main()
{
  std::ofstream logfile("output");
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TestPointValueHistory<2> test;
  test.run();
}
