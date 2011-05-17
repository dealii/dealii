//----------------------------  gerold_2.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2009 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  gerold_2.cc  ---------------------------


// apply SparsityTools::reorder_Cuthill_McKee to the cell connectivity
// graph for the mesh used in gerold_2. apparently the mesh consists
// of two or more non-connected parts, and the reordering algorithm
// trips over that

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <fstream>
#include <iomanip>

#include <deal.II/base/logstream.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria_boundary_lib.h>


template <int dim>
class LaplaceProblem
{
  public:
    void run ();
  private:
    Triangulation<dim>   triangulation;
};


template <int dim>
void LaplaceProblem<dim>::run ()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);

  std::ifstream input_file("gerold_1.inp");
  grid_in.read_ucd(input_file);

  SparsityPattern cell_connectivity;
  GridTools::get_face_connectivity_of_cells (triangulation,
					     cell_connectivity);
  std::vector<unsigned int> permutation(triangulation.n_active_cells());
  SparsityTools::reorder_Cuthill_McKee (cell_connectivity, permutation);

  for (unsigned int i=0; i<permutation.size(); ++i)
    deallog << permutation[i] << std::endl;
}


int main ()
{
  std::ofstream logfile("gerold_2/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      LaplaceProblem<3> laplace_problem_3d;
      laplace_problem_3d.run ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  return 0;
}
