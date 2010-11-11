// nedelec_2.cc,v 1.4 2003/04/09 15:49:55 wolf Exp
// (c) Wolfgang Bangerth
//
// Some more tests with the Nedelec element, basically using
// multi-element triangulations and check whether interface conditions
// are ok

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_tools.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_values.h>
#include <fe/mapping_cartesian.h>
#include <numerics/data_out.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2


void
plot_diff(const Vector<double>& v1, const Vector<double>& v2, const Tensor<1,2>& n)
{
  AssertDimension(v1.size(), 2);
  AssertDimension(v2.size(), 2);

  double p1 = v1(0)*n[1] - v1(1)*n[0];
  double p2 = v2(0)*n[1] - v2(1)*n[0];

  deallog << " tangential  diff " << p1-p2 << " (" << p1 << '-' << p2 << ')' << std::endl;
}


void
plot_diff(const Vector<double>& v1, const Vector<double>& v2, const Tensor<1,3>& n)
{
  AssertDimension(v1.size(), 3);
  AssertDimension(v2.size(), 3);
  
  Tensor<1,3> p1, p2;
  for (unsigned int d=0;d<3;++d)
    {
      const unsigned int d1 = (d+1)%3;
      const unsigned int d2 = (d+2)%3;
      
      p1[d] = v1(d1)*n[d2] - v1(d2)*n[d1];
      p2[d] = v2(d1)*n[d2] - v2(d2)*n[d1];
    }
  
  deallog << " tangential  diff " << p1-p2 << " (" << p1 << '-' << p2 << ')' << std::endl;
}


template <int dim>
void
plot (const Triangulation<dim> &tr, const unsigned int p)
{
  deallog << dim << "d, "
          << tr.n_active_cells() << " CELLS" << std::endl;

				   // Create the finite element space
  FE_Nedelec<dim> fe_ned (p);  
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe_ned);
  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof, cm);
  cm.close ();
  
                                   // generate some numbers for the
                                   // degrees of freedom on this mesh
  Vector<double> values (dof.n_dofs());
  for (unsigned int i=0; i<values.size(); ++i)
    values(i) = i;
  cm.distribute (values);
  
                                   // then make sure that hanging node
                                   // constraints are taken care of

                                   // now take these values, and print
                                   // the values of this so defined
                                   // function on each cell and on
                                   // each quadrature point
  QTrapez<dim-1> quadrature;
  std::vector<Vector<double> > shape_values1 (quadrature.size(),
                                             Vector<double>(dim));
  std::vector<Vector<double> > shape_values2 (quadrature.size(),
                                             Vector<double>(dim));
  
  MappingCartesian<dim> mapping;
  FEFaceValues<dim> feface(mapping, fe_ned, quadrature,
			   update_values|update_q_points|update_normal_vectors);
  FEFaceValues<dim> feneighbor(mapping, fe_ned, quadrature,
			       update_values|update_q_points|update_normal_vectors);
  FESubfaceValues<dim> fesubface(mapping, fe_ned, quadrature,
				 update_values|update_q_points|update_normal_vectors);
  
  
  for (typename DoFHandler<dim>::active_cell_iterator
         c = dof.begin_active();
       c!=dof.end(); ++c)
    for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
      {
	deallog << "cell " << c << " face " << face;
	if (c->at_boundary(face))
	  {
	    deallog << " boundary" << std::endl;
	    continue;
	  }

	if (c->neighbor(face)->has_children())
	  {
	    deallog << " neighbor " << c->neighbor(face) << " refined" << std::endl;
	    continue;
	  }
	
	deallog << " neighbor " << c->neighbor(face) << std::endl;
	
	feface.reinit(c, face);
	feface.get_function_values (values, shape_values1);
	
	if (c->neighbor(face)->level() < c->level())
	  {
	    std::pair<unsigned int, unsigned int> neighbor_face
	      = c->neighbor_of_coarser_neighbor(face);
	    
	    fesubface.reinit(c->neighbor(face), neighbor_face.first, neighbor_face.second);
	    fesubface.get_function_values (values, shape_values2);
	    for (unsigned int k=0;k<feface.n_quadrature_points;++k)
	      {
		    deallog << feface.quadrature_point(k);
		    plot_diff(shape_values1[k], shape_values2[k], feface.normal_vector(k));
	      }
	  }
	else
	  {
	    feneighbor.reinit(c->neighbor(face), c->neighbor_of_neighbor(face));
	    feneighbor.get_function_values (values, shape_values2);
	    for (unsigned int k=0;k<feface.n_quadrature_points;++k)
	      {
		deallog << feface.quadrature_point(k);
		plot_diff(shape_values1[k], shape_values2[k], feface.normal_vector(k));
	      }
	  }
      }
				   // Uncomment this to plot solutions
//   if (dim==2 && p==1)
//     {
//       std::ofstream gnuplot_output("test.gpl");
//       DataOut<dim> data_out;
//       data_out.attach_dof_handler (dof);
//       data_out.add_data_vector (values, "u");
//       data_out.build_patches (mapping, 4);
//       data_out.write_gnuplot(gnuplot_output);
//     }
}



template<int dim>
inline void
check (const unsigned int p)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

                                   // first everything on a
                                   // once-refined grid
  tr.refine_global (1);
//   plot (tr, p);

                                   // then same with one cell refined
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  plot (tr, p);
}


int
main()
{
  std::ofstream logfile ("nedelec_2/output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << "Degree 0:" << std::endl;
  check<2> (0);
  check<3> (0);
  deallog << "Degree 1:" << std::endl;
  check<2> (1);
  check<3> (1);
  
  return 0;
}



