
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


// check whether we can read a 3d grid. this test used to fail until late
// 2003, when the necessary infrastructure was created

#include "../tests.h"
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>


#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

#include <grid/grid_in.h>
#include <grid/tria_boundary_lib.h>

#include <base/timer.h>


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};




template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		dof_handler (triangulation)
{};


template <int dim>
void LaplaceProblem<dim>::run () 
{
  deallog << "Solving problem in " << dim << " space dimensions." << std::endl;
  
 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);


  std::ifstream input_file("gerold_1.inp");
 

  deallog << "read ucd data file" << std::endl;


  grid_in.read_ucd(input_file);
  deallog << "ucd data file readin exe" << std::endl;
};


int main () 
{
  std::ofstream logfile("gerold_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
				
  try
    {
      LaplaceProblem<3> laplace_problem_3d;
      laplace_problem_3d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
  
  return 0;
};
