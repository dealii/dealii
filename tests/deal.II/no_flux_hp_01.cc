// check the creation of no-flux boundary conditions for a finite
// element that consists of only a single set of vector components
// (i.e. it has dim components)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>


template<int dim>
void test (const Triangulation<dim>& tr,
		      const hp::FECollection<dim>& fe)
{
  hp::DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      deallog << "FE=" << fe[0].get_name()
	      << ", case=" << i
	      << std::endl;

      std::set<unsigned char> boundary_ids;
      for (unsigned int j=0; j<=i; ++j)
	boundary_ids.insert (j);
      
      ConstraintMatrix cm;
      VectorTools::compute_no_normal_flux_constraints (dof, 0, boundary_ids, cm);

      cm.print (deallog.get_file_stream ());
    }
}


template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    tr.begin_active()->face(i)->set_boundary_indicator (i);
  
  tr.refine_global(2);

  for (unsigned int degree=1; degree<4; ++degree)
    {
      hp::FECollection<dim> fe (FESystem<dim> (FE_Q<dim>(degree), dim));
      test(tr, fe);
    }
}


int main()
{
  std::ofstream logfile ("no_flux_hp_01/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
