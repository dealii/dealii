//----------------------------  point_value_history_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors and Michael Rapson
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point_value_history_01.cc  ---------------------------



// ************************************************
// A test program for the PointValueHistory class
// Currently this only tests a finite element system
// with 3 components, on a hyper cube and with Vectors.
// The testing is done in two dimensions.
// ************************************************




#include <base/quadrature_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_renumbering.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>

#include <lac/vector.h>

#include <numerics/point_value_history.h>

#include <fstream>

using namespace dealii;



template <int dim>
class TestPointValueHistory
{
public:
    TestPointValueHistory();
    void run();

private:
    Triangulation <dim> triangulation;
    FESystem<dim> finite_element;
    DoFHandler <dim > dof_handler;
    PointValueHistory <dim> test_copy;
};




template <int dim>
TestPointValueHistory<dim>::TestPointValueHistory() :
finite_element(FE_Q<dim > (1 + 1), dim, FE_Q<dim > (1), 1),

dof_handler(triangulation)
{ }




template <int dim>
void TestPointValueHistory<dim>::run()
{
    // Make a triangulation
    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(2); // refine 2 times to make 5 nodes per side

    // make a DOF handler, a model solution filled with ones and a flow vector
    FESystem<dim> finite_element(FE_Q<dim > (1 + 1), dim, FE_Q<dim > (1), 1);
    DoFHandler <dim > dof_handler(triangulation);
    dof_handler.distribute_dofs(finite_element);
    DoFRenumbering::Cuthill_McKee(dof_handler);

    // renumber for components so that same dof indices are used for BlockVectors and normal Vectors
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1; // component dim = pressure component!
    DoFRenumbering::component_wise(dof_handler, block_component);

    // Vector
    Vector <double> solution, post_processed, poles;
    solution.reinit(dof_handler.n_dofs());
    post_processed.reinit(dof_handler.n_dofs());
    poles.reinit(dof_handler.n_dofs());

    // set up a simple linear discrete time system so that time plots vary
    // over the mesh but can be easily checked. The basic idea is to have each
    // component of the fe_system to depend on a specific dimension (i.e component 0
    // depends on dim 0, etc. % dim handles the case where there are more components
    // than dimensions. The code breaks down at the edges of the mesh and this is
    // not corrected for.
    {
        QGauss<dim> quadrature_formula(2);
        FEValues<dim> fe_values(finite_element, quadrature_formula, update_values | update_quadrature_points); // just need local_dof_indices and quadrature_points

        std::vector<unsigned int> local_dof_indices(finite_element.dofs_per_cell);
        std::vector<Point<dim> > dof_locations(finite_element.dofs_per_cell);
        Vector<double> cell_pole(finite_element.dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell, endc;
        cell = dof_handler.begin_active();
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values.reinit(cell); // need to get local_dof_indices
            cell->get_dof_indices(local_dof_indices);
            dof_locations = fe_values.get_quadrature_points();
            cell_pole = 0;
            for (unsigned int dof = 0; dof != finite_element.dofs_per_cell; dof++)
            {
                unsigned int dof_component = finite_element.system_to_component_index(dof).first;

                for (unsigned int q_point = 0; q_point < quadrature_formula.size(); ++q_point)
                {
                    cell_pole(dof) += (fe_values.shape_value(dof, q_point) * dof_locations [q_point] (dof_component % dim));
                }
                solution(local_dof_indices [dof]) = 1; // start all solutions at 1
                poles(local_dof_indices[dof]) -= cell_pole(dof);

                if (dof_component == dim) // components start numbering at 0
                    poles(local_dof_indices[dof]) = -0.1; // dim+1th component is not handled well by the code above
            }

        } // loop over all cells
        poles.add(1.0); // slow down the pole settling time
        post_processed.add(3.0); // set to starting value.
    }

    // Setup monitor node to print variation over time
    unsigned int n_inputs = 1;
    PointValueHistory<dim> node_monitor (dof_handler, n_inputs);
    PointValueHistory<dim> no_dof_handler (n_inputs);

    // check that the assignment operator is valid
    test_copy = node_monitor;
    test_copy.add_point(Point<2>::Point(1, 0.2));
    test_copy.add_field_name("Solution");
    std::vector < std::vector <Point <dim> > > selected_locations;
    test_copy.get_points(selected_locations);
    test_copy.mark_locations();
    test_copy.close();
    test_copy.start_new_dataset(0.1);
    test_copy.evaluate_field("Solution", solution);
    std::vector <double> input_value(n_inputs, 1);
    test_copy.push_back_independent(input_value);
    test_copy.write_gnuplot("point_value_history_01/Test_Copy");
    test_copy.status(deallog.get_file_stream());
    test_copy.clear ();
    // end of assignment operator check

    {
        node_monitor.add_field_name("Solution");
        node_monitor.add_field_name("Post Processed Vector"); // not sensitive to spaces

        // two alternatives here, adding a point at a time or a vector of points
        // 2d points
        std::vector <Point < 2 > > point_vector(5, Point < 2 > ::Point());
        point_vector[0] = Point < 2 > ::Point(0, 0); // some of these points will hit a node, others won't
        point_vector[1] = Point < 2 > ::Point(0.25, 0);
        point_vector[2] = Point < 2 > ::Point(0.25, 0.45);
        point_vector[3] = Point < 2 > ::Point(0.45, 0.45);
        point_vector[4] = Point < 2 > ::Point(0.8, 0.8);

        node_monitor.add_points(point_vector);
        node_monitor.add_point(Point<2>::Point(1, 0.2)); // add a single point

        // MonitorNode requires that the instance is 'closed' before any data is added
        // this ensures that points are not added once time starts.
        node_monitor.close();
        no_dof_handler.close(); // closing still required!

        std::vector < std::vector <Point <dim> > > selected_locations;
        node_monitor.get_points(selected_locations);
	// write output to a file
        Vector<double> node_locations = node_monitor.mark_locations();
	// write output to a file
    }

    double delta_t = 0.000001;
    double t_max = 0.00001;

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

        solution.scale(poles); // decaying exponentials of varying time constants
        post_processed = solution;
        post_processed.add(2.0); // simple post processing, giving it a dc offset
    }
    node_monitor.write_gnuplot("point_value_history_01/node");
    no_dof_handler.write_gnuplot("point_value_history_01/no_dof");

    node_monitor.status (deallog.get_file_stream());
    no_dof_handler.status (deallog.get_file_stream());

    deallog << "Starting data files" << std::endl;

				     // copy all the data into deallog and
				     // delete those files
    const std::string filenames[]
      = { "point_value_history_01/node_00.gpl",
	  "point_value_history_01/node_01.gpl",
	  "point_value_history_01/node_02.gpl",
	  "point_value_history_01/node_03.gpl",
	  "point_value_history_01/node_04.gpl",
	  "point_value_history_01/node_05.gpl",
	  "point_value_history_01/node_indep.gpl",
	  "point_value_history_01/Test_Copy_00.gpl",
	  "point_value_history_01/Test_Copy_indep.gpl",
	  "point_value_history_01/no_dof_indep.gpl" };

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




int main()
{
    std::ofstream logfile("point_value_history_01/output");
    logfile << std::setprecision(2);
    deallog << std::setprecision(2);
    deallog.attach(logfile);
    deallog.depth_console(0);
    deallog.threshold_double(1.e-10);

    TestPointValueHistory<2> test;
    test.run();
}
