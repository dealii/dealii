// $Id$
// (c) Guido Kanschat
//
// Compute support points

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

template <int dim>
inline void
check_support (const FiniteElement<dim>& finel, const char* name)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (finel);

  std::cout << name << '<' << dim << '>' << " cell support points" << std::endl;
  
  const std::vector<Point<dim> > &cell_points = finel.get_unit_support_points ();

  for (unsigned int k=0;k<cell_points.size();++k)
    std::cout << std::setprecision(3) << cell_points[k] << std::endl;
  
  const std::vector<Point<dim-1> > &face_points = finel.get_unit_face_support_points ();
  const std::vector<double> dummy_weights (face_points.size());
  
  Quadrature<dim-1> q(face_points, dummy_weights);
  QProjector<dim> qp(q, false);
  
  unsigned int j=0;
  for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
    {
      std::cout << name << '<' << dim << '>' << " face " << i << " support points" << std::endl;
        
      for (unsigned int k=0;k<face_points.size();++k)
	std::cout << std::setprecision(3) << qp.point(j++)
	     << std::endl;
    }
}

template <int dim>
inline void
check_matrices (FiniteElement<dim>& fe, const char* name)
{
  std::cout << name << '<' << dim << '>' << " constraint " << std::endl;
  fe.constraints().print_formatted (std::cout, 7, false, 10, "~");

  for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
    {
      std::cout << name << '<' << dim << '>' << " restriction " << i << std::endl;
      fe.restrict(i).print_formatted (std::cout, 3, false, 6, "~");
      std::cout << name << '<' << dim << '>' << " embedding " << i << std::endl;
      fe.prolongate(i).print_formatted (std::cout, 3, false, 6, "~");
    }
}


#define CHECK_S(EL,deg,dim)   { FE_ ## EL<dim> EL(deg); check_support(EL, #EL #deg); }
#define CHECK_M(EL,deg,dim)   { FE_ ## EL<dim> EL(deg); check_matrices(EL, #EL #deg); }
#define CHECK_ALL(EL,deg,dim) { FE_ ## EL<dim> EL(deg); check_support(EL, #EL #deg); check_matrices(EL,#EL #deg); }

int
main()
{
  CHECK_M(DGQ,0,2);
  CHECK_M(DGQ,1,2);
  CHECK_M(DGQ,2,2);
  CHECK_M(DGQ,3,2);
  CHECK_M(DGQ,4,2);
  CHECK_ALL(Q,1,2);
  CHECK_ALL(Q,2,2);
  CHECK_ALL(Q,3,2);
  CHECK_ALL(Q,4,2);
  CHECK_M(DGQ,0,3);
  CHECK_M(DGQ,1,3);
  CHECK_M(DGQ,2,3);
  CHECK_M(DGQ,3,3);
  CHECK_M(DGQ,4,3);
  CHECK_ALL(Q,1,3);
  CHECK_ALL(Q,2,3);
//  CHECK_ALL(Q,3,2);
//  CHECK_ALL(Q,4,2);
  return 0;
}
