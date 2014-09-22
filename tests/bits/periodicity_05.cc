//----------------------------  periodicity_05.cc  ---------------------------
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2010, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  periodicity_05.cc  ---------------------------


// test for bug #82
// (http://code.google.com/p/dealii/issues/detail?id=82) which
// demonstrates that under some circumstances we create cycles in
// constraints

#include "../tests.h"
#include <iomanip>
#include <fstream>


#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/iterative_inverse.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/matrix_out.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

using namespace dealii;

class Deal2PeriodicBug
{

public:
  Deal2PeriodicBug();
  void run();
private:
  void makeGrid();
  void make_periodicity_constraints();
  void setup_system();

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;
  ConstraintMatrix constraints;
};

Deal2PeriodicBug::Deal2PeriodicBug()
:  fe(2), dof_handler(triangulation)
{}


void Deal2PeriodicBug::run()
{
  makeGrid();
  setup_system();
}

void Deal2PeriodicBug::make_periodicity_constraints()
{
  std::vector<bool> mask(1);
  mask[0] = true;
  ComponentMask cmask(mask);
  // Here we use the DoFTools function to place periodic constraints in the constraints matrix
  // We have set the boundary index for the left face of the boundary to 0 and the right face
  // index to 2 (Hence the 2nd and 3rd arguments) The direction indicator is to specify that
  // the DOFs are to be matched in the y-direction (0 in the 4th argument since they are allowed
  // to differ in the x-direction, it compares each of the coordinates not specified by the direction
  // integer which here is 0) Since we wanted to only specify the director and electric components
  // for periodicity we use a component mask as well
  DoFTools::make_periodicity_constraints(dof_handler, 2, 0, 0, constraints, cmask);
}

void Deal2PeriodicBug::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  deallog << "Making Constraint Matrix..." << std::endl;
  make_periodicity_constraints();
  deallog << "Constraint Matrix Complete" << std::endl;

  constraints.print(deallog.get_file_stream());

  constraints.close();
}

void Deal2PeriodicBug::makeGrid()
{
  deallog<< "Constructing the grid..." <<std::endl;
  const Point<2> p1(0,0), p2(1,1);
  GridGenerator::hyper_rectangle(triangulation,p1,p2);
  triangulation.begin_active()->face(2)->set_boundary_indicator(1);
  triangulation.begin_active()->face(3)->set_boundary_indicator(1);
  triangulation.begin_active()->face(0)->set_boundary_indicator(0);
  triangulation.begin_active()->face(1)->set_boundary_indicator(2);
  triangulation.refine_global(1);

  Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
  (++(++cell))->set_refine_flag();
  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  deallog<< "Number of active cells: " << triangulation.n_active_cells() << std::endl;
  GridOut grid_out;
  grid_out.write_eps (triangulation, deallog.get_file_stream());
  deallog<< "Grid construction complete..." <<std::endl;
}





int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Deal2PeriodicBug Deal2Bug;
  Deal2Bug.run();
  return 0;
}
