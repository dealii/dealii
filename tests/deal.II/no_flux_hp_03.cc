// check the creation of no-flux boundary conditions for a finite
// element that consists of only a single set of vector components
// (i.e. it has dim components)
//
// like no_flux_hp_01 but check on a hyper_sphere geometry

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>


template<int dim>
void test (const Triangulation<dim>& tr,
		      const hp::FECollection<dim>& fe)
{
  hp::DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe[0].get_name()
	  << std::endl;

  std::set<unsigned char> boundary_ids;
  boundary_ids.insert (0);
      
  ConstraintMatrix cm;
  VectorTools::compute_no_normal_flux_constraints (dof, 0, boundary_ids, cm);

  cm.print (deallog.get_file_stream ());
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);
  
  tr.refine_global(2);

  for (unsigned int degree=1; degree<4; ++degree)
    {
      hp::FECollection<dim> fe (FESystem<dim> (FE_Q<dim>(degree), dim));
      test(tr, fe);
    }
}


int main()
{
  std::ofstream logfile ("no_flux_hp_03/output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
