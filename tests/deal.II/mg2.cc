// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <numerics/multigrid.h>
#include <numerics/mg_smoother.h>
#include <grid/tria.h>
#include <grid/mg_dof.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>

main()
{
  Triangulation<2> tr;
  MGDoFHandler<2> dof(&tr);
  FELinear<2> fe;
  
  GridGenerator::hyper_cube(tr,-1.,1.);
  tr.refine_global(3);
  dof.distribute_dofs(fe);
}

