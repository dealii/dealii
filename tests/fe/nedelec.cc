// $Id$
// (c) Wolfgang Bangerth
//
// Show the shape functions of the Nedelec element

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2

char fname[50];

void
transform_grid (Triangulation<2> &tria,
		const unsigned int  transform)
{
  switch (transform)
    {
      case 0:
					     // first round: take
					     // original grid
	    break;
      case 1:
					     // second round: rotate
					     // triangulation
	    GridTools::rotate (3.14159265358/2, tria);
	    break;
      case 2:
					     // third round: inflate
					     // by a factor of 2
	    GridTools::scale (2, tria);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
};

	    


template<int dim>
inline void
plot_shape_functions()
{
  FE_Nedelec<dim> fe_ned(1);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

				   // check the following with a
				   // number of transformed
				   // triangulations
  for (unsigned int transform=0; transform<3; ++transform)
    {
      transform_grid (tr, transform);

      DoFHandler<dim> dof(tr);
      typename DoFHandler<dim>::cell_iterator c = dof.begin();
      dof.distribute_dofs(fe_ned);
      
      QTrapez<1> q_trapez;
      const unsigned int div=2;
      QIterated<dim> q(q_trapez, div);
      FEValues<dim> fe(fe_ned, q, update_values);
      fe.reinit(c);
      
      unsigned int k=0;
      for (unsigned int mz=0;mz<=((dim>2) ? div : 0) ;++mz)
	for (unsigned int my=0;my<=((dim>1) ? div : 0) ;++my)
	  for (unsigned int mx=0;mx<=div;++mx)
	    {
	      deallog << "q_point(" << k << ")=" << q.point(k)
		      << std::endl;
	      
	      for (unsigned int i=0;i<fe_ned.dofs_per_cell;++i)
		{
		  deallog << "  shape_function(" << i << ")=";
		  for (unsigned int c=0; c<fe.get_fe().n_components(); ++c)
		    deallog << " " << fe.shape_value_component(i,k,c);
		  deallog << std::endl;
		};
	      k++;
	    }
    }
};


int
main()
{
  std::ofstream logfile ("nedelec.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  plot_shape_functions<2>();
//  plot_shape_functions<3>();
  
  return 0;
}



