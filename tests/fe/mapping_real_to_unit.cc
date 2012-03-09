
#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q1.h>


template<int dim, int spacedim>
void test_real_to_unit_cell()
{
  
 Triangulation<dim, spacedim>   triangulation;
  GridGenerator::hyper_cube (triangulation);


  const unsigned int n_points = 5;
  
  
  std::vector< Point<dim> > unit_points(Utilities::fixed_power<dim>(n_points));


  switch (dim)
    {

      case 1:
	    for (unsigned int x=0; x<n_points;++x)
	      unit_points[x][0] = double(x)/double(n_points);
	    break;
	    
      case 2:
	    for (unsigned int x=0; x<n_points;++x)
	      for (unsigned int y=0; y<n_points;++y)
		{
		  unit_points[y * n_points + x][0] = double(x)/double(n_points);
		  unit_points[y * n_points + x][1] = double(y)/double(n_points);
		}
	    break;

      case 3:
	    for (unsigned int x=0; x<n_points;++x)
	      for (unsigned int y=0; y<n_points;++y)
		for (unsigned int z=0; z<n_points;++z)
		{
		  unit_points[z * n_points + y * n_points + x][0] = double(x)/double(n_points);
		  unit_points[z * n_points + y * n_points + x][1] = double(y)/double(n_points);
		  unit_points[z * n_points + y * n_points + x][2] = double(z)/double(n_points);
		}
	    break;
    }
  
    
  MappingQ1< dim, spacedim > map;

  
  typename Triangulation<dim, spacedim >::active_cell_iterator cell = triangulation.begin_active();

				   //Move a vertex a little bit
  const unsigned int n_dx = 5;
  const double dx = 0.4/n_dx;
  Point<spacedim> direction;
  for (unsigned int j=0; j<spacedim;++j)
    direction[j]=dx;
  
  for (unsigned int j=0; j<n_dx;++j)
    {
      cell->vertex(0) = double(j)*direction;
      deallog << "Vertex displacement: " << j*dx <<  std::endl <<  std::endl;
      for (unsigned int i=0; i<Utilities::fixed_power<dim>(n_points);++i)
	{
	  
	  Point<spacedim> p = map.transform_unit_to_real_cell(cell,unit_points[i]);
	  Point<dim> p_unit = map.transform_real_to_unit_cell(cell,p);
	  // deallog << "Orig   :" << unit_points[i] << std::endl;
	  // deallog << "p      :" << p << std::endl;
	  // deallog << "p_unit :" << p_unit <<  std::endl;
	  deallog << "Distance: " << unit_points[i].distance(p_unit) <<  std::endl <<  std::endl;
	}
      deallog << std::endl <<  std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("mapping_real_to_unit/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  test_real_to_unit_cell<1,1>();
  test_real_to_unit_cell<2,2>();
  test_real_to_unit_cell<3,3>();
  
  test_real_to_unit_cell<1,2>();
				   // test_real_to_unit_cell<1,3>();
  test_real_to_unit_cell<2,3>();
  
  return 0;
}



