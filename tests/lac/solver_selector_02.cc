//----------------------------------------------------------------------
//    $Id: solver_selector_01.cc 21777 2010-08-29 17:36:27Z bangerth $
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Test a bug in the SolverSelector when using a custom SolverControl. At one
// point the SolverControl got "sliced".

#include "../tests.h"
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_selector.h>

#include <fstream>

class MySolverControl:
   public SolverControl
{
  public:
    virtual State check(const unsigned int step,const double check_value)
      {
	deallog << "MySolverControl " << step << std::endl;
	return SolverControl::check(step, 0);
      }
};




template <class MATRIX, class VECTOR>
void
check(const MATRIX& A, const VECTOR& f)
{
  std::vector<std::string> names;
  names.push_back("cg");
  names.push_back("bicgstab");
  names.push_back("gmres");
  names.push_back("fgmres");

  MySolverControl mycont;
  SolverSelector<VECTOR> solver;
  PreconditionSSOR<SparseMatrix<double> > pre;
  pre.initialize(A);
  
  VECTOR u;
  u.reinit(f);
  
  std::vector<std::string>::const_iterator name;

  solver.set_control(mycont);  
  for (name = names.begin(); name != names.end();++name)
    {
      solver.select(*name);
      u = 0.;
      solver.solve(A, u, f, pre);
    }

				   //test deprecated constructor too:
  GrowingVectorMemory<VECTOR> mem;
  SolverSelector<VECTOR> solver2("gmres", mycont, mem);
  u = 0.;
  solver2.solve(A, u, f, pre);
}


int main()
{
  std::ofstream logfile("solver_selector_02/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  unsigned int size=37;
  unsigned int dim = (size-1)*(size-1);
  
  deallog << "Size " << size << " Unknowns " << dim << std::endl;
  
				   // Make matrix
  FDMatrix testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  testproblem.five_point(A);
  Vector<double>  f(dim);
  f = 1.;
  
  check(A, f);
}
