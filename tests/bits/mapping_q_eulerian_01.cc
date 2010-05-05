//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// test that DataOut and MappingQEulerian agree on how to output a
// displaced mesh. This was broken between r20158 and r21072

#include "../tests.h"

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <lac/full_matrix.h>
#include <lac/identity_matrix.h>
#include <base/quadrature_lib.h>
#include <base/thread_management.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/filtered_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/sparse_matrix.h>
#include <grid/grid_generator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_tools.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <lac/constraint_matrix.h>
#include <fe/mapping_q1.h>
#include <fe/mapping_q_eulerian.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>
#include <base/multithread_info.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>

using namespace dealii;


template <int dim>
class Displacement : public Function<dim>
{
  public:
    Displacement() :
		    Function<dim>(dim)
      {}

    double value (const Point<dim> &p,
		  const unsigned int component) const
      {
	return p[component];
      }

    void vector_value (const Point<dim> &p,
		       Vector<double> &v) const
      {
	for (unsigned int i=0; i<dim; ++i)
	  v(i) = p[i];
      }
};


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FESystem<dim> fe(FE_Q<dim>(1),dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> displacements (dof_handler.n_dofs());

  VectorTools::interpolate (dof_handler,
			    Displacement<dim>(),
			    displacements);

  MappingQEulerian<dim> euler(2, displacements, dof_handler);
				   // now the actual test
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  std::vector<std::string> names (dim, "displacement");
  data_out.add_data_vector(displacements,names);

  				   // output with all cells curved
  data_out.build_patches(euler,1,DataOut<dim>::curved_inner_cells);
  data_out.write_gnuplot(deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile("mapping_q_eulerian_01/output");
  deallog.attach(logfile);
  deallog << std::setprecision (4);
  logfile << std::setprecision (4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}


