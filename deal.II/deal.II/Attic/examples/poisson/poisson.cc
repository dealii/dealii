/* $Id$ */


#include <grid/dof.h>
#include <grid/tria.h>
#include <fe/fe_lib.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include "../../deal/dsmatrix.h>
#include <grid/dof_constraints.h>
#include <basic/data_io.h>
#include <numerics/problem_base.h>
#include <numerics/problem_assembler.h>

#include <fstream.h>
#include <cmath>
extern "C" {
#  include <stdlib.h>
}

extern TriaActiveIterator<1,CellAccessor<1> > x;
extern TriaActiveIterator<2,CellAccessor<2> > y;


void main () {
  Triangulation<1> tria;
  DoFHandler<1>    dof(tria);
  FELinear<1>      fe;
  ProblemBase<1>   problem(&tria, &dof);
  Quadrature<1>    q;
  
  tria.create_hypercube();
  dof.distribute_dofs (fe);
  problem.assemble (q);
};
